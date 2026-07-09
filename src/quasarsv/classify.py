"""Event classification — physiological IG/TR vs somatic translocation.

Set by `classify_events(calls)` after annotation. The classifications are:

* ``IG_intra``           — both ends in the same IG/TR locus.
                           V(D)J recombination or class-switch recombination.
                           **Physiological.**
* ``IG_IG``              — between two different IG/TR loci. Usually mapping
                           noise; occasionally legitimate (e.g. trans-rearr).
* ``IG_driver_canonical``— one end in an IG/TR locus, the other in an annotated
                           driver, AND the pair is in the canonical
                           translocation table (`data/known_partners.tsv`).
                           **Somatic, clinical.**
* ``IG_driver_novel``    — IG/TR ↔ annotated driver but NOT in the canonical
                           table — review-worthy putative driver event.
* ``IG_intergenic``      — IG/TR ↔ unannotated region. Often artefact, rarely
                           novel.
* ``driver_driver``      — two annotated drivers, different genes
                           (e.g. NPM1-ALK).  **Somatic, clinical.**
* ``driver_intra``       — within one driver (aSHM or short-range duplicate).
* ``driver_intergenic``  — driver ↔ unannotated region.
* ``intergenic``         — neither end annotated.

Helper `is_somatic_clinical(call)` returns True for the two clinically
actionable classes.
"""
from __future__ import annotations

from .annotate import load_builtin_loci, load_builtin_partners, _is_known_pair
from .model import FusionCall


IG_TR_LOCI = {
    "IGH", "IGK", "IGL", "IGH_Emu", "IGH_3RR",
    "TRA", "TRB", "TRG", "TRD",
}


def _is_driver(loci_by_gene: dict[str, str], gene: str) -> bool:
    role = loci_by_gene.get(gene, "")
    return role == "driver"


def _is_ig(gene: str) -> bool:
    return gene in IG_TR_LOCI


def classify_events(calls: list[FusionCall]) -> None:
    """Annotate every FusionCall with its event_class — mutates in place."""
    loci = load_builtin_loci()
    role = {g.gene: g.role for g in loci}
    partners = load_builtin_partners()

    for c in calls:
        ga, gb = c.gene_a, c.gene_b
        a_ig, b_ig = _is_ig(ga), _is_ig(gb)
        a_drv = _is_driver(role, ga)
        b_drv = _is_driver(role, gb)

        if a_ig and b_ig:
            if ga == gb:
                c.event_class = "IG_intra"
            else:
                c.event_class = "IG_IG"
        elif a_ig or b_ig:
            other = gb if a_ig else ga
            if _is_driver(role, other):
                if _is_known_pair(partners, ga, gb):
                    c.event_class = "IG_driver_canonical"
                else:
                    c.event_class = "IG_driver_novel"
            elif other:
                c.event_class = "IG_intergenic"
            else:
                c.event_class = "IG_intergenic"
        elif a_drv and b_drv:
            if ga == gb:
                c.event_class = "driver_intra"
            else:
                c.event_class = "driver_driver"
        elif a_drv or b_drv:
            c.event_class = "driver_intergenic"
        else:
            c.event_class = "intergenic"


SOMATIC_CLINICAL = {"IG_driver_canonical", "driver_driver"}
SOMATIC_WORTH_REVIEW = {"IG_driver_novel", "driver_intergenic"}
PHYSIOLOGICAL = {"IG_intra", "IG_IG"}


def is_somatic_clinical(c: FusionCall) -> bool:
    return c.event_class in SOMATIC_CLINICAL


def is_physiological(c: FusionCall) -> bool:
    return c.event_class in PHYSIOLOGICAL


def demote_nonclinical_t1(
    calls: list[FusionCall],
    target_tier: str = "T3",
) -> int:
    """Demote T1 calls that are physiological (V(D)J / class-switch),
    intra-gene (aSHM short-range duplicates), or flagged ``recurrent_artefact``.

    Single-caller forge_scan T1s reach this stage via the
    ``single_caller_very_strong`` rule; those that fall into one of these
    non-clinical classes pollute the T1 list and degrade precision against
    a canonical-translocation truth set without any clinical benefit.

    Returns the number of calls demoted.
    """
    demoted = 0
    for c in calls:
        if c.tier != "T1":
            continue
        if c.known_partner:
            # Always preserve known canonical partner pairs at their assigned tier.
            continue
        is_nonclinical_intra = c.event_class in ("IG_intra", "IG_IG", "driver_intra")
        is_recurrent_artefact = "recurrent_artefact" in c.qc_flags
        # Multi-caller PASS evidence overrides demotion — multi-caller agreement on
        # an intra-IG/driver hit is unusual and worth surfacing.
        if c.n_callers >= 2 and c.any_pass:
            continue
        if is_nonclinical_intra or is_recurrent_artefact:
            c.tier = target_tier
            flag = "demoted_nonclinical_t1"
            if flag not in c.qc_flags:
                c.qc_flags.append(flag)
            demoted += 1
    return demoted


def _has_two_genes(c: FusionCall) -> bool:
    return bool(c.gene_a) and bool(c.gene_b) and c.gene_a != c.gene_b


def gate_t1_actionable(calls: list[FusionCall]) -> int:
    """Cap T1 (clinically actionable) to a *resolved* gene-pair claim.

    A T1 call must be either a two-gene KNOWN canonical partner pair, or
    supported by ≥2 independent callers. A single-caller call whose partner
    side is unannotated (a lone driver breakpoint such as BCL6→intergenic, or a
    lone IG breakpoint) is review-worthy but not actionable, so it is capped to
    T2. This preserves recall (T2 still counts) while keeping the T1 list to
    breakpoints a reporting scientist can act on.

    Returns the number of calls capped.
    """
    capped = 0
    for c in calls:
        if c.tier != "T1":
            continue
        canonical = _has_two_genes(c) and c.known_partner
        multicaller = c.n_callers >= 2
        if not (canonical or multicaller):
            c.tier = "T2"
            if "t1_capped_unresolved" not in c.qc_flags:
                c.qc_flags.append("t1_capped_unresolved")
            capped += 1
    return capped
