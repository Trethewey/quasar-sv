"""Post-merge QC flags — detect recurrent-position artefacts and other red flags."""
from __future__ import annotations

import csv
import re
from collections import Counter
from importlib.resources import files as pkg_files
from pathlib import Path

from .model import FusionCall

# Non-primary GRCh38 contigs: decoys, unplaced/unlocalised scaffolds, ALT
# haplotypes, HLA, EBV. A breakpoint whose partner lands here is essentially
# always a mapping artefact, never a real clinical translocation partner.
_DECOY_CONTIG_RE = re.compile(
    r"(_decoy|_random|_alt|^chrUn|^Un_|KI270|GL000|JTFH|KN707|HLA-|EBV)", re.I)


def _load_artefact_loci() -> list[tuple[str, int, int, str]]:
    p = Path(str(pkg_files("quasarsv").joinpath("data/artefact_loci.tsv")))
    out: list[tuple[str, int, int, str]] = []
    if not p.exists():
        return out
    with open(p, encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            chrom = row["chrom"]
            chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
            out.append((chrom, int(row["start"]), int(row["end"]), row.get("notes", "")))
    return out


def flag_builtin_artefact_loci(calls: list[FusionCall]) -> None:
    """Flag any breakpoint that lands in the curated artefact-loci table.

    These are well-characterised mapping noise hotspots that produce false
    partner calls regardless of caller; they should always be suspect.
    """
    art = _load_artefact_loci()
    if not art:
        return
    for c in calls:
        for chrom, start, end, _ in art:
            if (c.chrom_a == chrom and start <= c.pos_a <= end) or \
               (c.chrom_b == chrom and start <= c.pos_b <= end):
                if "builtin_artefact_locus" not in c.qc_flags:
                    c.qc_flags.append("builtin_artefact_locus")
                # downgrade — these are NEVER real calls
                if c.tier == "T1":
                    c.tier = "T3"
                elif c.tier == "T2":
                    c.tier = "T3"
                break


def flag_recurrent_position_artefacts(
    calls: list[FusionCall],
    window: int = 500,
    min_partners: int = 10,
    min_distinct_chroms: int = 3,
    auto_downgrade: bool = False,
) -> None:
    """Flag breakpoints where many distinct partners on many chromosomes hit
    the same ±`window` bp.

    Real biology like promiscuous BCL6 has many partners but usually clustered
    on a few chromosomes; library/mapping artefacts produce mates spanning the
    whole genome. We therefore require BOTH a high partner count AND mates on
    ≥`min_distinct_chroms` chromosomes before flagging.

    Advisory by default — set `auto_downgrade=True` to drop T1→T2 on flag.
    """
    partner_positions: dict[tuple[str, int], set[tuple[str, int]]] = {}
    for c in calls:
        ka = (c.chrom_a, c.pos_a // window)
        kb = (c.chrom_b, c.pos_b // window)
        partner_positions.setdefault(ka, set()).add((c.chrom_b, c.pos_b // window))
        partner_positions.setdefault(kb, set()).add((c.chrom_a, c.pos_a // window))

    hot = set()
    for k, partners in partner_positions.items():
        if len(partners) >= min_partners:
            chroms = {p[0] for p in partners}
            if len(chroms) >= min_distinct_chroms:
                hot.add(k)

    for c in calls:
        ka = (c.chrom_a, c.pos_a // window)
        kb = (c.chrom_b, c.pos_b // window)
        if ka in hot or kb in hot:
            if "recurrent_artefact" not in c.qc_flags:
                c.qc_flags.append("recurrent_artefact")
            if auto_downgrade:
                if c.tier == "T1":
                    c.tier = "T2"
                elif c.tier == "T2":
                    c.tier = "T3"


# Known aSHM (activation-induced cytidine deaminase) target genes in B-cell
# lymphomas — short-range clustered breaks here are usually biology, not noise.
ASHM_TARGETS = {
    "BCL6", "BCL7A", "BTG1", "BTG2", "MYC", "PAX5", "PIM1", "CXCR4",
    "IRF4", "RHOH", "ST6GAL1", "SOCS1", "REL", "CIITA",
}


def flag_short_range_intrachr(
    calls: list[FusionCall],
    min_span: int = 200,
    exempt_genes: set[str] | None = None,
) -> None:
    """Mark intra-chr breakpoints within `min_span` bp as 'short_range'.

    Skips calls in known aSHM target genes (these are real biology — clustered
    short-range breakpoints reflect AID-driven hypermutation, not artefact).
    """
    exempt = ASHM_TARGETS if exempt_genes is None else exempt_genes
    for c in calls:
        if c.gene_a in exempt or c.gene_b in exempt:
            continue
        if c.chrom_a == c.chrom_b and abs(c.pos_b - c.pos_a) < min_span:
            c.qc_flags.append("short_range")


def cohort_recurrence(calls: list[FusionCall], window: int = 1000) -> Counter:
    """Recurrence across samples by (chrom_a, pos_a//window, chrom_b, pos_b//window)."""
    key = lambda c: (c.chrom_a, c.pos_a // window, c.chrom_b, c.pos_b // window)
    seen: dict[tuple, set[str]] = {}
    for c in calls:
        seen.setdefault(key(c), set()).add(c.sample)
    return Counter({k: len(s) for k, s in seen.items() if len(s) > 1})


def _is_decoy_partner(c: FusionCall) -> bool:
    return bool(_DECOY_CONTIG_RE.search(c.chrom_a)) or bool(_DECOY_CONTIG_RE.search(c.chrom_b))


def flag_decoy_partner(calls: list[FusionCall]) -> int:
    """Downgrade to T3 any call with a breakpoint on a decoy / unplaced / ALT /
    HLA / EBV contig. These are mapping artefacts, never real clinical partners.

    A genuine two-gene canonical known-partner pair is exempt (defensive — a
    canonical pair should never have a decoy endpoint, but do not silently bury
    one if annotation somehow produced it).

    Returns the number of calls downgraded.
    """
    n = 0
    for c in calls:
        if c.tier == "T3":
            continue
        if c.known_partner and c.gene_a and c.gene_b:
            continue
        if _is_decoy_partner(c):
            c.tier = "T3"
            if "decoy_partner" not in c.qc_flags:
                c.qc_flags.append("decoy_partner")
            n += 1
    return n


# Event classes that carry no clinical signal on their own — physiological
# recombination and intra-locus duplicates. Downgrading these out of the review
# tier keeps the T1/T2 list dominated by candidate somatic events.
_PHYSIOLOGICAL_CLASSES = {"IG_intra", "IG_IG", "driver_intra"}
# Classes we must NEVER auto-downgrade on an artefact flag alone — a novel
# driver rearrangement is exactly what a lymphoma-specific caller should surface
# (see audit: demoting these on recurrent_artefact buries real novel drivers).
_CLINICAL_CLASSES = {"IG_driver_canonical", "IG_driver_novel", "driver_driver", "driver_intergenic"}


def demote_physiological_noise(calls: list[FusionCall]) -> int:
    """Downgrade physiological and recurrent-artefact NON-clinical calls out of
    the T1/T2 review tier to T3.

    Removed from the review tier:
      * physiological classes (IG_intra / IG_IG / driver_intra), and
      * ``recurrent_artefact``-flagged calls whose class is IG_intergenic or
        intergenic (the IG-locus-to-everywhere mapping noise).

    Preserved (never downgraded here):
      * known canonical partner pairs, and
      * clinically relevant classes (canonical / novel driver / driver-driver /
        driver-intergenic) — so novel drivers are not buried by an artefact flag.

    Returns the number of calls downgraded.
    """
    n = 0
    for c in calls:
        if c.tier not in ("T1", "T2"):
            continue
        if c.known_partner or c.event_class in _CLINICAL_CLASSES:
            continue
        physiological = c.event_class in _PHYSIOLOGICAL_CLASSES
        recurrent_noise = ("recurrent_artefact" in c.qc_flags
                           and c.event_class in ("IG_intergenic", "intergenic"))
        if physiological or recurrent_noise:
            c.tier = "T3"
            if "demoted_review_noise" not in c.qc_flags:
                c.qc_flags.append("demoted_review_noise")
            n += 1
    return n


def apply_default_qc(
    calls: list[FusionCall],
    sample_lineage: dict[str, str] | None = None,
    lineage_default: str = "B",
) -> list[FusionCall]:
    """One-shot: built-in artefact mask + recurrent-position + short-range flags
    + IG-switch-region rescue synthesis + canonical-partner tier promotion.

    Parameters
    ----------
    calls
        FusionCall list (mutated in place).
    sample_lineage
        Optional ``{sample_id: "B" | "T" | "any"}`` overriding the rescue's
        lineage prior per sample.
    lineage_default
        Lineage assumption used when ``sample_lineage`` is empty / missing the
        sample. Defaults to ``"B"`` (B-cell lymphoma — the package's primary
        use-case).
    """
    flag_builtin_artefact_loci(calls)
    flag_recurrent_position_artefacts(calls)
    flag_short_range_intrachr(calls)
    # Rescue runs AFTER masking so it can see which calls were artefact-flagged.
    from .rescue import RescueConfig, rescue_ig_driver_pairs
    rescue_ig_driver_pairs(
        calls,
        cfg=RescueConfig(lineage=lineage_default),
        sample_lineage=sample_lineage,
    )
    # Final pass: canonical-partner promotion (after rescue has added gene_a/gene_b)
    from .promote import promote_known_partners
    promote_known_partners(calls)
    # Event classification — physiological vs somatic
    from .classify import classify_events, demote_nonclinical_t1
    classify_events(calls)
    # Final precision pass: demote single-caller T1 calls that are
    # physiological (IG_intra / IG_IG), intra-driver (aSHM duplicate), or
    # flagged ``recurrent_artefact``. Multi-caller PASS and known-canonical
    # pairs are preserved.
    demote_nonclinical_t1(calls)
    # Review-tier hygiene: remove mapping-artefact and physiological noise from
    # T1/T2 so the review list stays dominated by candidate somatic events.
    flag_decoy_partner(calls)
    demote_physiological_noise(calls)
    # Actionability gate: T1 is reserved for resolved canonical gene pairs or
    # multi-caller support; lone-gene single-caller breakpoints cap at T2.
    from .classify import gate_t1_actionable
    gate_t1_actionable(calls)
    # NOTE: empirical-LLR promotion (src/quasarsv/llr_score.py) is shipped
    # but NOT wired in by default — testing it on the cohort dropped F1 from
    # 0.84 to 0.32 because the heuristic Phred-per-read calibration is too
    # generous for quasarsv's evidence distribution (it promoted 53 noisy
    # calls). The module + tests remain available for downstream callers who
    # want to opt in via a custom apply_default_qc(... + promote_by_llr).
    return calls
