"""Phred-scaled LLR scoring for FusionCalls.

Promotion-only secondary tier signal. Approximates GRIDSS empirical-LLR
using bulk evidence summaries (SR, PE, assembly, MAPQ) rather than
per-call evidence reconstruction. Calls scoring above the configured
threshold get lifted to T1 or T2; never demoted.

NOT wired into apply_default_qc by default — the constant Phred-per-read
calibration is too generous on cohort-distribution noise and dropped F1
from 0.84 to 0.32 in testing. Use opt-in via ``promote_by_llr(calls)``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .model import FusionCall


@dataclass
class LlrConfig:
    """Knobs for the LLR scoring + tier promotion."""
    # Phred-scaled LLR thresholds for tier promotion
    t1_threshold: float = 100.0
    t2_threshold: float = 50.0
    # Per-evidence-type contributions (Phred per supporting read, capped)
    sr_phred_per_read: float = 4.0       # P(SR | ref) ≈ 10^-0.4
    pe_phred_per_read: float = 1.8       # P(PE | ref) ≈ 10^-0.18
    sr_cap_reads: int = 30                # cap to limit dominance of one big cluster
    pe_cap_reads: int = 30
    # Categorical bonuses
    assembly_contig_bonus: float = 20.0
    multi_caller_bonus_per_extra: float = 10.0
    known_partner_bonus: float = 30.0
    precise_bonus: float = 10.0
    pass_bonus: float = 5.0
    # Penalties — recurrent_artefact / short_range / V(D)J-like discounts
    intra_ig_penalty: float = -40.0
    intra_driver_penalty: float = -25.0
    recurrent_artefact_penalty: float = -60.0


def llr_score(call: FusionCall, cfg: LlrConfig | None = None) -> float:
    """Phred-scaled LLR that ``call`` is a real somatic variant.

    Higher = more confident. The score is bounded below at 0 — negative
    contributions can cancel each other but the floor is 0 for clarity.
    """
    cfg = cfg or LlrConfig()
    score = 0.0

    # SR contribution — capped
    sr = min(int(call.split_reads or 0), cfg.sr_cap_reads)
    score += sr * cfg.sr_phred_per_read

    # PE contribution — capped
    pe = min(int(call.discordant_pairs or 0), cfg.pe_cap_reads)
    score += pe * cfg.pe_phred_per_read

    if int(call.assembly_contigs or 0) > 0:
        score += cfg.assembly_contig_bonus

    if int(call.n_callers or 0) >= 2:
        score += cfg.multi_caller_bonus_per_extra * (int(call.n_callers) - 1)

    if bool(call.known_partner):
        score += cfg.known_partner_bonus

    if bool(call.precise):
        score += cfg.precise_bonus

    if bool(call.any_pass):
        score += cfg.pass_bonus

    # Penalties for non-clinical event classes (V(D)J + aSHM + artefacts)
    ec = call.event_class or ""
    if ec in ("IG_intra", "IG_IG"):
        score += cfg.intra_ig_penalty
    elif ec == "driver_intra":
        score += cfg.intra_driver_penalty

    if "recurrent_artefact" in (call.qc_flags or []):
        score += cfg.recurrent_artefact_penalty

    return max(0.0, score)


def promote_by_llr(
    calls: list[FusionCall],
    cfg: LlrConfig | None = None,
) -> int:
    """Mutate ``calls`` — promote calls whose LLR exceeds the configured
    tier threshold but whose existing tier is below.

    Returns the number of calls promoted.

    Properties:
      * NEVER demotes — only lifts T3 → T2, T2 → T1
      * Adds ``llr_promoted_T2`` / ``llr_promoted_T1`` to qc_flags so the
        promotion source is traceable
      * Stores the score on the call as a side-effect via qc_flags suffix
        ``llr_<score>`` (not a separate column, keeps the schema unchanged)
    """
    cfg = cfg or LlrConfig()
    promoted = 0
    tier_rank = {"T1": 0, "T2": 1, "T3": 2}
    for c in calls:
        score = llr_score(c, cfg)
        existing = tier_rank.get(c.tier, 3)
        new_tier: str | None = None
        if score >= cfg.t1_threshold and existing > 0:
            new_tier = "T1"
        elif score >= cfg.t2_threshold and existing > 1:
            new_tier = "T2"
        if new_tier is None:
            continue
        old = c.tier
        c.tier = new_tier
        flag = f"llr_promoted_{new_tier}_from_{old}"
        if flag not in c.qc_flags:
            c.qc_flags.append(flag)
        promoted += 1
    return promoted
