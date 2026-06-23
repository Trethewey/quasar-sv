"""Tier promotion for known canonical partner pairs.

Default tier rules are evidence-strict; they correctly avoid false positives
but they also bury T3 calls whose gene pair matches a *canonical lymphoma
translocation* (IGH-BCL2, MYC-IGH, IGH-BCL6, IGH-CCND1, NPM1-ALK …).

Clinical reading benefits from a separate, gene-pair-aware promotion: any
fusion call whose annotated gene pair is in `data/known_partners.tsv` and
which has any non-trivial evidence (`SR + PE >= 5` by default) is bumped to
at least T2. Pure split-read or precise breakpoints earn T1.
"""
from __future__ import annotations

from dataclasses import dataclass

from .model import FusionCall


@dataclass
class PromotionConfig:
    min_evidence_any: int = 5    # SR + PE total to qualify for promotion
    promote_to_t1_min_sr: int = 5
    promote_to_t1_min_pe: int = 10


def promote_known_partners(
    calls: list[FusionCall],
    cfg: PromotionConfig | None = None,
) -> int:
    """Mutate `calls` in place: promote known canonical partners based on evidence.

    Returns number of calls promoted.
    """
    cfg = cfg or PromotionConfig()
    promoted = 0
    tier_rank = {"T1": 0, "T2": 1, "T3": 2}
    for c in calls:
        if not c.known_partner:
            continue
        ev_sum = c.split_reads + c.discordant_pairs
        if ev_sum < cfg.min_evidence_any:
            continue
        # Strong: promote to T1 (e.g., IGH-MYC with SR>=5 OR PE>=10)
        if c.split_reads >= cfg.promote_to_t1_min_sr or c.discordant_pairs >= cfg.promote_to_t1_min_pe:
            target = "T1"
        else:
            target = "T2"
        if tier_rank.get(c.tier, 3) > tier_rank.get(target, 3):
            c.tier = target
            if "promoted_known_partner" not in c.qc_flags:
                c.qc_flags.append("promoted_known_partner")
            promoted += 1
    return promoted
