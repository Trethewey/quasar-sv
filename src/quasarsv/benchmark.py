"""Truth-set loading and gene-pair-level precision/recall scoring.

Truth-set TSV schema (see ``data/cohort_truth.tsv``):

    sample_id      matches the FusionCall.sample value
    cell_line      human-readable
    cohort         PMBL / DLBCL_ABC / DLBCL_GCB / DH_DLBCL / ATLL / ...
    truth_class    confirmed | likely | none_expected | driver_focal
    gene_a, gene_b expected partner gene pair (order-insensitive); '-' if N/A
    cytoband       e.g. t(3;14)(q27;q32)
    source         citation / DOI / cell-line bank accession
    notes          free text

``truth_class = none_expected`` is a negative control: any T1 call in that
sample counts as a false positive. Multiple truth rows may share a sample_id
(e.g. double-hit DLBCL with two canonical translocations).
"""
from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .model import FusionCall


@dataclass
class TruthEntry:
    sample_id: str
    cell_line: str = ""
    cohort: str = ""
    truth_class: str = "confirmed"          # confirmed | likely | none_expected | driver_focal
    gene_a: str = ""
    gene_b: str = ""
    cytoband: str = ""
    source: str = ""
    notes: str = ""

    @property
    def pair(self) -> frozenset[str]:
        """Order-insensitive gene pair (empty strings filtered)."""
        return frozenset(g for g in (self.gene_a, self.gene_b) if g and g != "-")

    @property
    def is_positive(self) -> bool:
        return self.truth_class in ("confirmed", "likely", "driver_focal")

    @property
    def is_negative(self) -> bool:
        return self.truth_class == "none_expected"


def load_truth_set(path: str | Path) -> list[TruthEntry]:
    """Read a truth-set TSV. Returns one TruthEntry per row."""
    p = Path(path)
    out: list[TruthEntry] = []
    if not p.exists():
        return out
    with open(p, encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            out.append(TruthEntry(
                sample_id=row.get("sample_id", "").strip(),
                cell_line=row.get("cell_line", "").strip(),
                cohort=row.get("cohort", "").strip(),
                truth_class=row.get("truth_class", "confirmed").strip() or "confirmed",
                gene_a=row.get("gene_a", "").strip(),
                gene_b=row.get("gene_b", "").strip(),
                cytoband=row.get("cytoband", "").strip(),
                source=row.get("source", "").strip(),
                notes=row.get("notes", "").strip(),
            ))
    return out


def _call_pair(c: FusionCall) -> frozenset[str]:
    return frozenset(g for g in (c.gene_a, c.gene_b) if g)


@dataclass
class SampleScore:
    sample: str
    expected_pairs: list[frozenset[str]] = field(default_factory=list)
    is_negative_control: bool = False
    # TP / FN per expected pair
    matched_pairs: list[frozenset[str]] = field(default_factory=list)
    missed_pairs: list[frozenset[str]] = field(default_factory=list)
    # Tier at which each matched pair surfaced (best-of)
    match_tier: dict[frozenset[str], str] = field(default_factory=dict)
    # All T1 calls in this sample (gene-pair frozensets)
    t1_call_pairs: list[frozenset[str]] = field(default_factory=list)
    # FP: T1 calls whose gene pair matches no truth entry (and the sample is
    # not specifically a negative control)
    fp_t1_count: int = 0


@dataclass
class BenchmarkScore:
    per_sample: dict[str, SampleScore]
    tool_name: str = "quasarsv"
    # Aggregates over samples with truth
    tp: int = 0
    fp: int = 0
    fn: int = 0

    @property
    def precision(self) -> float:
        denom = self.tp + self.fp
        return self.tp / denom if denom else 0.0

    @property
    def recall(self) -> float:
        denom = self.tp + self.fn
        return self.tp / denom if denom else 0.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if (p + r) else 0.0


def score_calls_against_truth(
    calls: list[FusionCall],
    truth: list[TruthEntry],
    tool_name: str = "quasarsv",
    accept_tiers: tuple[str, ...] = ("T1", "T2"),
    ambiguous_alts_count_as_fp: bool = False,
    restrict_to_scored_samples: bool = True,
    relax_canonical_ig_partner: bool = False,
    canonical_ig_set: frozenset[str] = frozenset(
        {"IGH", "IGK", "IGL", "IGH_Emu", "IGH_3RR",
         "TRA", "TRB", "TRG", "TRD"}),
) -> BenchmarkScore:
    """Score `calls` against `truth`.

    A truth entry matches a call when their gene-pair sets are equal (order
    insensitive). A sample's TP for that truth = the best (lowest-tier) call
    that matches; FN = no matching call at an accepted tier.

    A sample with ``truth_class = none_expected`` (negative control) has
    zero TPs / FNs by definition; every T1 call in that sample contributes
    to FP.

    For samples WITH expected truths, any T1 call whose gene pair is not in
    the sample's expected set increments FP. ``ig_partner_ambiguous`` calls
    can optionally be excluded from FP counting via
    ``ambiguous_alts_count_as_fp=False`` (the default — the rescue's
    intentional ambiguous alternatives shouldn't penalise precision).

    Parameters
    ----------
    restrict_to_scored_samples
        Skip truth entries whose sample appears in no call. Useful when the
        truth set anticipates future scans so "unscanned" doesn't conflate
        with "missed". Default True.
    relax_canonical_ig_partner
        Treat a driver-IG canonical pair as a match for a truth driver-IG
        canonical pair as long as both IG sides are in
        ``canonical_ig_set``. Lets BCL6-IGL count as a match for truth
        BCL6-IGH (both represent the same driver rearrangement to *an* IG
        locus). Default False (strict gene-pair match).
    """
    truth_by_sample: dict[str, list[TruthEntry]] = defaultdict(list)
    for t in truth:
        truth_by_sample[t.sample_id].append(t)

    calls_by_sample: dict[str, list[FusionCall]] = defaultdict(list)
    for c in calls:
        calls_by_sample[c.sample].append(c)

    per_sample: dict[str, SampleScore] = {}
    tier_rank = {"T1": 0, "T2": 1, "T3": 2}
    accept_rank = {tier_rank[t] for t in accept_tiers}

    def _match(call_pair: frozenset[str], expected_pair: frozenset[str]) -> bool:
        if call_pair == expected_pair:
            return True
        if relax_canonical_ig_partner:
            # Driver = non-IG gene; require driver equality + both sides have
            # an IG locus partner. This collapses BCL6-IGL ≡ BCL6-IGH for
            # truth-set matching when both are canonical IG partners.
            c_drv = call_pair - canonical_ig_set
            e_drv = expected_pair - canonical_ig_set
            c_ig = call_pair & canonical_ig_set
            e_ig = expected_pair & canonical_ig_set
            if c_drv and c_drv == e_drv and c_ig and e_ig:
                return True
            # Driver-only call vs driver-IG truth: count as match when the
            # driver is the same and the call has the driver as its sole
            # annotated gene (typical when the IG side of an artefact-routed
            # call has no gene annotation). Honest: we identified the driver
            # rearrangement, just couldn't resolve the partner.
            if (e_drv and c_drv == e_drv and not c_ig
                    and len(call_pair) == 1):
                return True
        return False

    for sample, sample_truth in truth_by_sample.items():
        sample_calls = calls_by_sample.get(sample, [])
        if restrict_to_scored_samples and not sample_calls:
            continue   # skip unscanned samples
        score = SampleScore(sample=sample)
        score.is_negative_control = all(t.is_negative for t in sample_truth)
        score.expected_pairs = [t.pair for t in sample_truth if t.is_positive]

        # T1 calls in this sample (gene-pair sets)
        for c in sample_calls:
            if c.tier == "T1":
                score.t1_call_pairs.append(_call_pair(c))

        # Match each expected pair against the call set
        for exp in score.expected_pairs:
            best_tier: str | None = None
            for c in sample_calls:
                if tier_rank.get(c.tier, 9) not in accept_rank:
                    continue
                if _match(_call_pair(c), exp):
                    if best_tier is None or tier_rank[c.tier] < tier_rank[best_tier]:
                        best_tier = c.tier
            if best_tier is not None:
                score.matched_pairs.append(exp)
                score.match_tier[exp] = best_tier
            else:
                score.missed_pairs.append(exp)

        # FP: T1 calls whose gene-pair matches no expected pair
        for c in sample_calls:
            if c.tier != "T1":
                continue
            if not ambiguous_alts_count_as_fp and "ig_partner_ambiguous" in c.qc_flags:
                continue
            cp = _call_pair(c)
            if not cp:
                continue
            if not any(_match(cp, exp) for exp in score.expected_pairs):
                score.fp_t1_count += 1

        per_sample[sample] = score

    bench = BenchmarkScore(per_sample=per_sample, tool_name=tool_name)
    for s in per_sample.values():
        bench.tp += len(s.matched_pairs)
        bench.fn += len(s.missed_pairs)
        bench.fp += s.fp_t1_count
    return bench


def write_benchmark_tsv(score: BenchmarkScore, path: str | Path) -> None:
    """Write per-sample TP/FP/FN + the aggregate P/R/F1 line."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "expected_pairs", "matched_pairs", "missed_pairs",
                    "match_tiers", "t1_calls", "fp_t1", "negative_control"])
        for sample, s in sorted(score.per_sample.items()):
            w.writerow([
                sample,
                ";".join("-".join(sorted(p)) for p in s.expected_pairs),
                ";".join("-".join(sorted(p)) for p in s.matched_pairs),
                ";".join("-".join(sorted(p)) for p in s.missed_pairs),
                ";".join(f"{'-'.join(sorted(p))}@{s.match_tier[p]}"
                         for p in s.matched_pairs),
                len(s.t1_call_pairs),
                s.fp_t1_count,
                "yes" if s.is_negative_control else "no",
            ])
        w.writerow([])
        w.writerow([f"# tool={score.tool_name}",
                    f"TP={score.tp}", f"FP={score.fp}", f"FN={score.fn}",
                    f"precision={score.precision:.4f}",
                    f"recall={score.recall:.4f}",
                    f"f1={score.f1:.4f}"])


def builtin_truth_path() -> Path:
    from importlib.resources import files as pkg_files
    return Path(str(pkg_files("quasarsv").joinpath("data/cohort_truth.tsv")))
