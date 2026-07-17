"""Evidence-level merger.

Clusters `BreakpointCall`s across callers into `FusionCandidate`s using
single-linkage on breakpoint pair position (within tolerance, matching strands).
Then assigns a tier based on independent evidence types rather than naive
caller-count.

Design rules
------------
* Two callers seeing the same split-read cluster do NOT count as two
  evidence types — see `FusionCandidate.evidence_summary`.
* Caller-count is recorded but NOT the tiering primitive.
* Tier rules are explicit and configurable via `TierThresholds`.
"""
from __future__ import annotations

from dataclasses import dataclass

from .model import BreakpointCall, FusionCandidate, FusionCall


@dataclass
class MergeConfig:
    """Breakpoint-clustering tolerances."""
    pos_tolerance: int = 250        # bp window for treating two breakpoints as the same
    same_strand_required: bool = True
    # Clustering algorithm:
    #   "bucket"  = fixed-grid bucket + neighbour scan (legacy; bucket-edge artefacts)
    #   "dbscan"  = density-based clustering: two points join if both pos_a AND pos_b
    #               are within pos_tolerance; single-linkage transitive closure.
    #               Eliminates bucket-boundary false splits — see
    #               docs/precision_techniques.md technique #7.
    clustering: str = "dbscan"
    dbscan_min_samples: int = 1     # 1 keeps singletons (current behaviour);
                                    # raise to 2 to require ≥2 supporting calls per cluster


@dataclass
class TierThresholds:
    """How calls map to tiers based on evidence.

    Evidence-type counts within a single caller are NOT treated as truly
    independent (they typically share underlying reads). T1 therefore requires
    BOTH multi-caller agreement AND multi-evidence-type support.

    T1 (high-confidence — clinically actionable):
        ≥2 callers supporting AND
        ≥1 caller with FILTER=PASS AND
        ≥3 split reads (max across callers) AND
        ≥2 independent evidence types AND
        breakpoint precise in ≥1 caller.
    T2 (moderate — review-worthy):
        ≥2 callers, OR
        single caller PASS with ≥5 split reads AND assembly contig.
    T3:
        anything else surviving sanity filters.
    """
    t1_min_callers: int = 2
    t1_require_pass: bool = True
    t1_min_split_reads: int = 3
    t1_min_evidence_types: int = 2
    t1_require_precise: bool = True

    t2_min_callers: int = 2
    t2_single_min_split_reads: int = 5
    t2_single_require_assembly: bool = True

    min_quality_for_any_tier: float = 0.0
    drop_singleton_imprecise: bool = True


def _key(call: BreakpointCall, tol: int) -> tuple:
    # Bucket position to tolerance so clusters within ±tol on a side become
    # candidate matches; final acceptance uses real distance.
    return (call.chrom_a, call.chrom_b, call.strand_a, call.strand_b,
            call.pos_a // tol, call.pos_b // tol)


def _close(a: BreakpointCall, b: BreakpointCall, tol: int, require_strand: bool) -> bool:
    if a.chrom_a != b.chrom_a or a.chrom_b != b.chrom_b:
        return False
    if require_strand and (a.strand_a != b.strand_a or a.strand_b != b.strand_b):
        return False
    return abs(a.pos_a - b.pos_a) <= tol and abs(a.pos_b - b.pos_b) <= tol


def _cluster_bucket(calls: list[BreakpointCall], cfg: MergeConfig) -> list[FusionCandidate]:
    """Legacy bucket-based single-linkage clustering. Retained for compat."""
    buckets: dict[tuple, list[BreakpointCall]] = {}
    for c in calls:
        k = _key(c, cfg.pos_tolerance)
        buckets.setdefault(k, []).append(c)

    used = [False] * len(calls)
    by_idx = {id(c): i for i, c in enumerate(calls)}
    candidates: list[FusionCandidate] = []

    for i, c in enumerate(calls):
        if used[i]:
            continue
        used[i] = True
        members = [c]
        chrom_a, chrom_b, sa, sb, ba, bb = _key(c, cfg.pos_tolerance)
        neigh = []
        for dba in (-1, 0, 1):
            for dbb in (-1, 0, 1):
                neigh.extend(buckets.get((chrom_a, chrom_b, sa, sb, ba + dba, bb + dbb), []))
        for o in neigh:
            j = by_idx[id(o)]
            if used[j] or j == i:
                continue
            if _close(c, o, cfg.pos_tolerance, cfg.same_strand_required):
                used[j] = True
                members.append(o)

        members.sort(key=lambda m: (m.pos_a, m.pos_b))
        rep = members[len(members) // 2]
        candidates.append(FusionCandidate(
            sample=c.sample,
            chrom_a=rep.chrom_a, pos_a=rep.pos_a, strand_a=rep.strand_a,
            chrom_b=rep.chrom_b, pos_b=rep.pos_b, strand_b=rep.strand_b,
            sv_type=rep.sv_type,
            evidences=[m.evidence for m in members],
            member_callers=sorted({m.evidence.caller for m in members}),
            member_record_ids=[m.record_id for m in members],
        ))
    return candidates


def _cluster_dbscan(calls: list[BreakpointCall], cfg: MergeConfig) -> list[FusionCandidate]:
    """Density-based clustering (DBSCAN-style with Chebyshev distance).

    Two calls are direct neighbours when:
        same chrom_a, chrom_b
        same strand_a, strand_b (if cfg.same_strand_required)
        |pos_a_i - pos_a_j| <= eps  AND  |pos_b_i - pos_b_j| <= eps

    Clusters = transitive closure of the neighbour relation. Eliminates the
    bucket-edge artefacts where a true cluster straddles two adjacent
    pos_tolerance buckets.

    Implementation: union-find over per-group calls. We still group by
    (chrom_a, chrom_b, [strand_a, strand_b]) first so the pairwise search
    is small. Inside a group, we sort by pos_a and only test pairs whose
    pos_a delta is ≤ eps (sliding window) — keeps complexity near O(n)
    on realistic data even though worst case is O(n²).
    """
    eps = cfg.pos_tolerance
    require_strand = cfg.same_strand_required

    # Group by (chrom_a, chrom_b [, strand_a, strand_b])
    groups: dict[tuple, list[int]] = {}
    for i, c in enumerate(calls):
        key = (c.chrom_a, c.chrom_b)
        if require_strand:
            key = key + (c.strand_a, c.strand_b)
        groups.setdefault(key, []).append(i)

    # Union-find
    parent = list(range(len(calls)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    for indices in groups.values():
        if len(indices) == 1:
            continue
        # Sort by pos_a; sliding window over pos_a within eps, then linear
        # filter on pos_b.
        indices.sort(key=lambda i: calls[i].pos_a)
        n = len(indices)
        left = 0
        for right in range(n):
            ci = calls[indices[right]]
            while left < right and ci.pos_a - calls[indices[left]].pos_a > eps:
                left += 1
            for j in range(left, right):
                cj = calls[indices[j]]
                if abs(ci.pos_b - cj.pos_b) <= eps:
                    union(indices[right], indices[j])

    # Aggregate by root
    cluster_members: dict[int, list[int]] = {}
    for i in range(len(calls)):
        root = find(i)
        cluster_members.setdefault(root, []).append(i)

    candidates: list[FusionCandidate] = []
    min_samples = max(1, cfg.dbscan_min_samples)
    for root, idxs in cluster_members.items():
        if len(idxs) < min_samples:
            continue
        members = [calls[i] for i in idxs]
        members.sort(key=lambda m: (m.pos_a, m.pos_b))
        rep = members[len(members) // 2]
        candidates.append(FusionCandidate(
            sample=rep.sample,
            chrom_a=rep.chrom_a, pos_a=rep.pos_a, strand_a=rep.strand_a,
            chrom_b=rep.chrom_b, pos_b=rep.pos_b, strand_b=rep.strand_b,
            sv_type=rep.sv_type,
            evidences=[m.evidence for m in members],
            member_callers=sorted({m.evidence.caller for m in members}),
            member_record_ids=[m.record_id for m in members],
        ))
    return candidates


def cluster_calls(calls: list[BreakpointCall], cfg: MergeConfig | None = None) -> list[FusionCandidate]:
    cfg = cfg or MergeConfig()
    if cfg.clustering == "bucket":
        return _cluster_bucket(calls, cfg)
    return _cluster_dbscan(calls, cfg)


def assign_tier(cand: FusionCandidate, t: TierThresholds | None = None) -> tuple[str, list[str]]:
    """Return (tier, qc_flags)."""
    t = t or TierThresholds()
    flags: list[str] = []
    summary = cand.evidence_summary()
    n_ev_types = cand.n_independent_evidence_types()
    n_callers = len(cand.member_callers)
    any_pass = cand.any_filter_pass()
    any_precise = cand.any_precise()
    max_sr = summary["split_read"]
    max_asm = summary["assembly_contig"]
    max_q = max((e.raw_qual for e in cand.evidences), default=0.0)

    if t.drop_singleton_imprecise and n_callers == 1 and not any_precise:
        flags.append("singleton_imprecise")

    if max_q < t.min_quality_for_any_tier:
        flags.append("low_quality")

    # T1 — multi-caller required
    if (n_callers >= t.t1_min_callers
            and (any_pass or not t.t1_require_pass)
            and max_sr >= t.t1_min_split_reads
            and n_ev_types >= t.t1_min_evidence_types
            and (any_precise or not t.t1_require_precise)
            and "singleton_imprecise" not in flags):
        return "T1", flags

    # Single-caller exceptional path: when one caller (typically quasar)
    # reports overwhelming SR+PE support with both evidence types, escalate.
    very_strong_single = (
        n_callers == 1 and any_pass and any_precise
        and max_sr >= 10 and summary["discordant_pair"] >= 10
        and n_ev_types >= 2
    )
    if very_strong_single:
        flags.append("single_caller_very_strong")
        return "T1", flags

    # T2 — multi-caller, or single caller with very strong PASS + assembly
    single_caller_strong = (
        n_callers == 1 and any_pass
        and max_sr >= t.t2_single_min_split_reads
        and (max_asm >= 1 or not t.t2_single_require_assembly)
    )
    if n_callers >= t.t2_min_callers or single_caller_strong:
        if n_callers == 1:
            flags.append("single_caller_strong")
        return "T2", flags

    # Single-caller medium support (still useful at known-partner sites).
    if n_callers == 1 and any_pass and max_sr >= 5:
        flags.append("single_caller_medium")
        return "T2", flags

    return "T3", flags


def candidates_to_calls(candidates: list[FusionCandidate], t: TierThresholds | None = None) -> list[FusionCall]:
    """Project clustered candidates into the final `FusionCall` shape (pre-annotation)."""
    out: list[FusionCall] = []
    for idx, cand in enumerate(candidates):
        tier, flags = assign_tier(cand, t)
        s = cand.evidence_summary()
        fusion_id = f"{cand.sample}__{cand.chrom_a}_{cand.pos_a}_{cand.strand_a}__{cand.chrom_b}_{cand.pos_b}_{cand.strand_b}__{idx:05d}"
        out.append(FusionCall(
            sample=cand.sample,
            fusion_id=fusion_id,
            chrom_a=cand.chrom_a, pos_a=cand.pos_a, strand_a=cand.strand_a,
            chrom_b=cand.chrom_b, pos_b=cand.pos_b, strand_b=cand.strand_b,
            sv_type=cand.sv_type,
            callers_supporting=cand.member_callers,
            n_callers=len(cand.member_callers),
            split_reads=s["split_read"],
            discordant_pairs=s["discordant_pair"],
            assembly_contigs=s["assembly_contig"],
            soft_clips=s["soft_clip"],
            n_evidence_types=cand.n_independent_evidence_types(),
            vaf=cand.max_vaf(),
            precise=cand.any_precise(),
            any_pass=cand.any_filter_pass(),
            raw_qual_max=max((e.raw_qual for e in cand.evidences), default=0.0),
            tier=tier,
            qc_flags=flags,
            member_record_ids=cand.member_record_ids,
        ))
    return out


def merge_caller_calls(
    per_caller: dict[str, list[BreakpointCall]],
    cfg: MergeConfig | None = None,
    tiers: TierThresholds | None = None,
) -> list[FusionCall]:
    """Top-level entry: takes a {caller: [BreakpointCall]} mapping, returns merged FusionCall list."""
    all_calls: list[BreakpointCall] = []
    for caller, calls in per_caller.items():
        all_calls.extend(calls)
    cands = cluster_calls(all_calls, cfg)
    return candidates_to_calls(cands, tiers)
