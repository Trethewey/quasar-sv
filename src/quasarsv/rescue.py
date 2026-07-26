"""Artefact-masked breakend handling.

This module used to *synthesise* driver-IG fusion calls: whenever a driver locus
and an IG locus both sent split reads to the GRCh38 poly-G attractor at
chr2:32,916, it paired them using a canonical-partner lookup plus a B-cell
lineage prior, and promoted canonical pairs to T1 with the artefact-side read
counts displayed as if they were junction support.

That inference was invalid and has been removed. Measurements on the WGS
validation cohort:

* **The artefact channel carries no signal.** Every locus sheds reads to
  chr2:32,916 at an essentially constant rate — ~200-280 split reads per 10k
  reads — whether or not it is rearranged. In Karpas-1106P the rate for BCL6
  (204/10k, rearrangement claimed) is *lower* than for TP53 (240/10k) and
  NOTCH1 (337/10k), neither of which is rearranged. "Both loci hit the same
  artefact" is therefore true of every pair of loci in the genome and says
  nothing about partnership.
* **The reads are not genomic.** The clipped segments feeding the attractor are
  Illumina adapter read-through (``AGATCGGAAGAGC``) and 2-colour-chemistry
  poly-G tails, not IG switch sequence. They are filtered at source now — see
  ``scanners.cram_scanner.is_noise_clip``.
* **The masking excuse does not hold.** Discordant mates do not depend on the
  junction sequence, so a poly-G attractor cannot hide them. MD903 shows
  BCL6-IGH at PE=16 by discordant pairs alone. Karpas-1106P and U2940, whose
  BCL6-IGH the rescue used to synthesise, have *zero* BCL6-IGH reads of either
  kind. Their translocation is absent from the data, not masked.

A real driver-IG translocation is recovered by the ordinary scanner path, from
reads that actually join the two loci. What cannot be measured is not reported:
there is deliberately no inference-based fallback here, because the channel it
would have to infer from is sequencing noise.

What remains is honest annotation: mark calls whose partner breakend falls in a
masked artefact region so that downstream tiering can refuse to treat an
unresolvable breakend as a resolved fusion.
"""
from __future__ import annotations

from dataclasses import dataclass

from .model import FusionCall


# Antigen-receptor locus sets. Retained because event classification and
# reporting still need to know which loci are IG/TR; no longer used to choose a
# partner. ``_allowed_igs_for_lineage`` is deliberately gone with the prior.
BCELL_IGS = {"IGH", "IGK", "IGL", "IGH_Emu", "IGH_3RR"}
TCELL_IGS = {"TRA", "TRB", "TRG", "TRD"}
IG_LOCI = BCELL_IGS | TCELL_IGS


@dataclass
class RescueConfig:
    """Configuration for artefact-masked breakend annotation.

    The partner-inference knobs (``min_artefact_sr_per_side``, ``ratio_keep``,
    ``emit_canonical_alternatives``, ``noncanonical_fanout_ratio``, the pair
    caps) are gone: they tuned a fabrication. No threshold on a signal-free
    channel can produce a valid partner assignment.

    ``lineage`` is gone too, and it is worth being explicit about why: the
    B-cell/T-cell prior existed solely to restrict which IG locus the rescue was
    allowed to NAME as a partner. With the naming removed, the prior has no
    consumer — nothing downstream reads it. Keeping the field would advertise a
    control that does nothing.
    """
    # Tier ceiling for a call whose partner breakend is unresolvable.
    unresolved_tier: str = "T3"


def flag_artefact_masked_breakends(
    calls: list[FusionCall],
    cfg: RescueConfig | None = None,
) -> list[FusionCall]:
    """Annotate calls whose partner breakend lands in a masked artefact region.

    Such a call localises one real breakend but cannot name what it joins to.
    It is marked ``partner_undetermined`` and capped at a review tier so it can
    never be presented as a resolved fusion. No partner is invented, and no
    synthetic call is created.

    Returns the same list, mutated in place.
    """
    cfg = cfg or RescueConfig()
    tier_rank = {"T1": 0, "T2": 1, "T3": 2}
    cap = tier_rank.get(cfg.unresolved_tier, 2)

    for c in calls:
        if "builtin_artefact_locus" not in c.qc_flags:
            continue
        # One breakend sits in a masked artefact: whatever this locus joins to,
        # the artefact side cannot identify it.
        if "partner_undetermined" not in c.qc_flags:
            c.qc_flags.append("partner_undetermined")
        if tier_rank.get(c.tier, 2) < cap:
            c.tier = cfg.unresolved_tier
    return calls


def rescue_ig_driver_pairs(*_args, **_kwargs):
    """Removed: fabricated IG partner assignments from a signal-free channel.

    See the module docstring. Use :func:`flag_artefact_masked_breakends` for
    honest annotation of unresolvable breakends; use the scanner for detection.
    """
    raise NotImplementedError(
        "rescue_ig_driver_pairs has been removed: it assigned IG partners from a "
        "canonical lookup plus a lineage prior, with no read linking driver to "
        "partner, and promoted the result to T1. The artefact channel it relied on "
        "is adapter/poly-G sequencing noise emitted uniformly by every locus. "
        "Use flag_artefact_masked_breakends() instead."
    )
