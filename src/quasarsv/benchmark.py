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
        """A scoreable positive needs a gene pair.

        A row naming no genes (``gene_a = '-'``) cannot be matched by any
        gene-pair call, so counting it as a positive would charge every tool a
        false negative it had no way to avoid. Such rows document a sample
        (e.g. a published SV burden with no named partners) without scoring it.
        """
        return (self.truth_class in ("confirmed", "likely", "driver_focal")
                and bool(self.pair))

    @property
    def is_negative(self) -> bool:
        return self.truth_class == "none_expected"

    @property
    def is_disputed(self) -> bool:
        """Quarantined: literature and reads disagree, the material's identity is
        contested, or the event is outside the panel's vocabulary.

        A disputed row is scored neither as a true positive nor a false negative,
        and a call matching it is not a false positive. Resolving such a row in
        whichever direction flattered the tool is exactly the error this
        benchmark exists to avoid — so it is excluded, and counted in the
        report so the exclusion stays visible.
        """
        return self.truth_class == "disputed"


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
    # Detection-vs-lookup split: matched pairs for which a direct read-level
    # junction exists in the CRAM, versus those the tool got right without one.
    detected_pairs: list[frozenset[str]] = field(default_factory=list)
    lookup_only_pairs: list[frozenset[str]] = field(default_factory=list)


@dataclass
class BenchmarkScore:
    per_sample: dict[str, SampleScore]
    tool_name: str = "quasarsv"
    # Aggregates over samples with truth
    tp: int = 0
    fp: int = 0
    fn: int = 0
    # Of `tp`, how many were backed by a direct read-level junction (detection)
    # versus how many the tool named without one (lookup / inference). Only
    # populated when a junction-support set is supplied.
    tp_detected: int = 0
    tp_lookup_only: int = 0
    # Visibility of what was NOT scored. Quarantining contested truth is
    # defensible; doing it silently is not — a reader must be able to see how
    # much of the cohort the headline actually covers.
    disputed_rows: int = 0          # quarantined truth rows encountered
    skipped_samples: list[str] = field(default_factory=list)   # unscoreable
    fp_exempted_disputed: int = 0   # calls matching a quarantined row

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
    fp_tiers: tuple[str, ...] | None = None,
    ambiguous_alts_count_as_fp: bool = True,
    restrict_to_scored_samples: bool = True,
    relax_canonical_ig_partner: bool = False,
    match_driver_only: bool = False,
    junction_support: set[tuple[str, frozenset[str]]] | None = None,
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
    the sample's expected set increments FP.

    Scoring honesty
    ---------------
    The headline metric must measure DETECTION, not agreement between a
    canonical-partner prior and a literature-derived truth set that encodes the
    same canonical translocations. Three former leniencies each made a call
    match without the tool having resolved the fusion, and all three are now
    off by default:

    * ``relax_canonical_ig_partner`` — let BCL6-IGL match truth BCL6-IGH.
    * ``match_driver_only`` — let a driver-only call match a driver-IG truth.
    * ``ambiguous_alts_count_as_fp=False`` — exempted a tool's own hedged
      alternates from precision. Emitting IGH, IGK *and* IGL for one driver and
      being scored only on the one that hits is not precision.

    Enable them only to report a clearly-labelled lenient figure ALONGSIDE the
    strict one, never as the headline.

    Parameters
    ----------
    fp_tiers
        Tiers at which a non-truth call counts as a false positive. Defaults to
        ``accept_tiers``, so a tool is charged for exactly the tiers it is
        credited for. Counting TPs at T1+T2 while charging FPs only at T1 gives
        free recall to whichever tool emits the most T2 — a structural reward
        for volume, not accuracy, worth ~12k calls to one tool here and ~500 to
        another. Pass ``("T1",)`` to report the clinically-actionable-only view
        alongside the review-list view.
    restrict_to_scored_samples
        Skip truth entries whose sample appears in no call. Useful when the
        truth set anticipates future scans so "unscanned" doesn't conflate
        with "missed". Default True.
    relax_canonical_ig_partner
        Treat a driver-IG canonical pair as a match for a truth driver-IG
        canonical pair as long as both IG sides are in ``canonical_ig_set``.
        Default False (strict gene-pair match).
    match_driver_only
        Treat a call annotating only the driver as a match for a driver-IG
        truth. Requires ``relax_canonical_ig_partner``. Default False.
    junction_support
        ``{(sample, frozenset({gene_a, gene_b}))}`` for which a direct
        read-level junction was independently established in the reads. When
        supplied, each TP is split into detected (junction present) versus
        lookup-only (named correctly without one). A tool scoring TPs with no
        junction support is reproducing a prior, not detecting an event.
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
    fp_tier_set = set(fp_tiers if fp_tiers is not None else accept_tiers)
    skipped: list[str] = []
    n_disputed = 0
    n_fp_exempt = 0

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
            # Driver-only call vs driver-IG truth. "We found the driver but
            # could not resolve the partner" is a PARTIAL result, not a
            # detected fusion; counting it as a match makes a tool that never
            # resolves a partner score identically to one that always does.
            # Off unless explicitly requested for a lenient side-metric.
            if (match_driver_only and e_drv and c_drv == e_drv and not c_ig
                    and len(call_pair) == 1):
                return True
        return False

    for sample, sample_truth in truth_by_sample.items():
        sample_calls = calls_by_sample.get(sample, [])
        if restrict_to_scored_samples and not sample_calls:
            continue   # skip unscanned samples
        n_disputed += sum(1 for t in sample_truth if t.is_disputed)
        # A sample with no positive and no negative truth is not scoreable at
        # all: entirely-disputed, or documented-only (a published SV burden with
        # no named gene pair). Without this it would present as an empty expected
        # set, which silently turns every call into a false positive — the
        # opposite of quarantine. Recorded so the exclusion is visible: these
        # samples leave FP accounting too, which is a real limitation of the
        # headline, not a neutral omission.
        if not any(t.is_positive or t.is_negative for t in sample_truth):
            skipped.append(sample)
            continue
        score = SampleScore(sample=sample)
        score.is_negative_control = all(t.is_negative for t in sample_truth
                                        if not t.is_disputed)
        score.expected_pairs = [t.pair for t in sample_truth if t.is_positive]
        disputed_pairs = [t.pair for t in sample_truth if t.is_disputed and t.pair]

        # T1 calls in this sample (gene-pair sets)
        for c in sample_calls:
            if c.tier == "T1":
                score.t1_call_pairs.append(_call_pair(c))
        disputed_sets = [d for d in disputed_pairs if d]

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
                if junction_support is not None:
                    if (sample, exp) in junction_support:
                        score.detected_pairs.append(exp)
                    else:
                        score.lookup_only_pairs.append(exp)
            else:
                score.missed_pairs.append(exp)

        # FP: calls at an FP-charged tier whose gene-pair matches no expected pair
        for c in sample_calls:
            if c.tier not in fp_tier_set:
                continue
            if not ambiguous_alts_count_as_fp and "ig_partner_ambiguous" in c.qc_flags:
                continue
            cp = _call_pair(c)
            if not cp:
                continue
            if any(_match(cp, exp) for exp in score.expected_pairs):
                continue
            # A call matching a quarantined (disputed) truth row is neither
            # credited nor charged — the evidence is genuinely unresolved.
            if any(cp == d or _match(cp, d) for d in disputed_sets):
                n_fp_exempt += 1
                continue
            score.fp_t1_count += 1

        per_sample[sample] = score

    bench = BenchmarkScore(per_sample=per_sample, tool_name=tool_name)
    for s in per_sample.values():
        bench.tp += len(s.matched_pairs)
        bench.fn += len(s.missed_pairs)
        bench.fp += s.fp_t1_count
        bench.tp_detected += len(s.detected_pairs)
        bench.tp_lookup_only += len(s.lookup_only_pairs)
    bench.disputed_rows = n_disputed
    bench.skipped_samples = sorted(skipped)
    bench.fp_exempted_disputed = n_fp_exempt
    return bench


def write_benchmark_tsv(score: BenchmarkScore, path: str | Path) -> None:
    """Write per-sample TP/FP/FN + the aggregate P/R/F1 line."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        # NB "fp_calls" counts calls at `fp_tiers` (default = accept_tiers, i.e.
        # T1+T2), not T1 only — the column was previously named fp_t1 while
        # counting both.
        w.writerow(["sample", "expected_pairs", "matched_pairs", "missed_pairs",
                    "match_tiers", "t1_calls", "fp_calls", "negative_control",
                    "detected_pairs", "lookup_only_pairs"])
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
                ";".join("-".join(sorted(p)) for p in s.detected_pairs),
                ";".join("-".join(sorted(p)) for p in s.lookup_only_pairs),
            ])
        w.writerow([])
        w.writerow([f"# tool={score.tool_name}",
                    f"TP={score.tp}", f"FP={score.fp}", f"FN={score.fn}",
                    f"precision={score.precision:.4f}",
                    f"recall={score.recall:.4f}",
                    f"f1={score.f1:.4f}",
                    f"tp_detected={score.tp_detected}",
                    f"tp_lookup_only={score.tp_lookup_only}"])
        # What the headline does NOT cover, stated rather than omitted.
        w.writerow([f"# excluded: disputed_truth_rows={score.disputed_rows}",
                    f"fp_exempted_by_quarantine={score.fp_exempted_disputed}",
                    f"unscoreable_samples={len(score.skipped_samples)}",
                    ";".join(score.skipped_samples)])


def load_junction_support(
    path: str | Path,
    min_sr: int = 0,
    min_pe: int = 0,
) -> set[tuple[str, frozenset[str]]]:
    """Load independently-established read-level junctions.

    Expects a TSV with ``sample``, ``gene_a``, ``gene_b``. If a ``supported``
    column is present ("yes"/"true"/"1"), it is authoritative and the count
    thresholds are ignored — the caller has already applied a background model.
    Otherwise a row qualifies when it clears ``min_sr`` OR ``min_pe``.

    Returns ``{(sample, frozenset({gene_a, gene_b}))}`` for use as
    ``score_calls_against_truth(junction_support=...)``.
    """
    p = Path(path)
    out: set[tuple[str, frozenset[str]]] = set()
    if not p.exists():
        return out
    with open(p, encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            sample = (row.get("sample") or "").strip()
            ga = (row.get("gene_a") or "").strip()
            gb = (row.get("gene_b") or "").strip()
            if not sample or not ga or not gb:
                continue
            if "supported" in (r.fieldnames or []):
                if (row.get("supported") or "").strip().lower() not in ("yes", "true", "1"):
                    continue
            else:
                sr = int(float(row.get("sr_junction") or 0))
                pe = int(float(row.get("pe_junction") or 0))
                if sr < min_sr and pe < min_pe:
                    continue
            out.add((sample, frozenset({ga, gb})))
    return out


def builtin_truth_path() -> Path:
    from importlib.resources import files as pkg_files
    return Path(str(pkg_files("quasarsv").joinpath("data/cohort_truth.tsv")))
