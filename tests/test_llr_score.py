"""Tests for the empirical-LLR scoring + tier promotion."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.llr_score import LlrConfig, llr_score, promote_by_llr
from quasarsv.model import FusionCall


def _fc(tier="T3", sr=0, pe=0, asm=0, n_callers=1,
        known_partner=False, precise=False, any_pass=False,
        event_class="", qc_flags=None):
    return FusionCall(
        sample="S1", fusion_id="f1",
        chrom_a="14", pos_a=106_000_000, strand_a="+",
        chrom_b="18", pos_b=63_200_000, strand_b="-",
        sv_type="BND",
        tier=tier, split_reads=sr, discordant_pairs=pe,
        assembly_contigs=asm, n_callers=n_callers,
        known_partner=known_partner, precise=precise, any_pass=any_pass,
        event_class=event_class, qc_flags=list(qc_flags or []),
    )


def test_llr_is_nonnegative():
    call = _fc()
    assert llr_score(call) >= 0


def test_llr_grows_with_evidence():
    a = _fc(sr=2, pe=2)
    b = _fc(sr=10, pe=10)
    assert llr_score(b) > llr_score(a)


def test_llr_caps_runaway_sr():
    """A single huge SR cluster shouldn't dominate the score."""
    huge = _fc(sr=10_000)
    bounded = _fc(sr=30)
    # cap = 30 → both should score about the same on the SR contribution
    assert llr_score(huge) == llr_score(bounded)


def test_llr_known_partner_bonus():
    no_kp = _fc(sr=5, pe=10)
    kp = _fc(sr=5, pe=10, known_partner=True)
    assert llr_score(kp) > llr_score(no_kp)


def test_llr_multi_caller_bonus():
    single = _fc(sr=5, pe=10, n_callers=1)
    multi = _fc(sr=5, pe=10, n_callers=4)
    assert llr_score(multi) > llr_score(single)


def test_llr_intra_ig_penalty():
    """V(D)J IG_intra calls get a steep penalty even with high evidence."""
    vdj = _fc(sr=20, pe=10, event_class="IG_intra")
    real = _fc(sr=20, pe=10, event_class="IG_driver_canonical")
    assert llr_score(vdj) < llr_score(real)


def test_llr_recurrent_artefact_penalty():
    art = _fc(sr=20, pe=10, qc_flags=["recurrent_artefact"])
    clean = _fc(sr=20, pe=10)
    assert llr_score(art) < llr_score(clean)


def test_promote_lifts_t3_to_t2():
    # SR=10×4=40 + PE=10×1.8=18 + precise=10 + any_pass=5 = 73, in [50,100) → T2
    call = _fc(tier="T3", sr=10, pe=10, precise=True, any_pass=True)
    n = promote_by_llr([call])
    assert n == 1
    assert call.tier == "T2"
    assert any("llr_promoted_T2_from_T3" in f for f in call.qc_flags)


def test_promote_lifts_to_t1_when_score_high():
    # Multi-caller + assembly + known_partner + lots of SR/PE → > 100
    call = _fc(tier="T3", sr=20, pe=20, asm=1, n_callers=3,
               known_partner=True, precise=True, any_pass=True)
    n = promote_by_llr([call])
    assert n == 1
    assert call.tier == "T1"


def test_promote_never_demotes():
    """Even if LLR is low, an existing higher tier is never lowered."""
    call = _fc(tier="T1", sr=1, pe=0)
    promote_by_llr([call])
    assert call.tier == "T1"


def test_promote_no_change_when_already_above_threshold():
    """T1 call with high LLR stays T1, not flagged again."""
    call = _fc(tier="T1", sr=20, pe=20, known_partner=True)
    n = promote_by_llr([call])
    assert n == 0
    assert call.tier == "T1"
    assert not any(f.startswith("llr_promoted") for f in call.qc_flags)


def test_promote_idempotent():
    """Calling twice doesn't double-flag."""
    call = _fc(tier="T3", sr=10, pe=10, known_partner=True, precise=True, any_pass=True)
    promote_by_llr([call])
    flags_after_first = list(call.qc_flags)
    promote_by_llr([call])
    assert call.qc_flags == flags_after_first


def test_llr_thresholds_configurable():
    """Custom thresholds change which calls are promoted."""
    call = _fc(tier="T3", sr=5, pe=5)
    # Default thresholds — likely won't promote
    base_n = promote_by_llr([_fc(tier="T3", sr=5, pe=5)])
    # Very low thresholds — will promote
    cfg = LlrConfig(t1_threshold=1.0, t2_threshold=0.5)
    n = promote_by_llr([call], cfg)
    assert n == 1
    assert call.tier in ("T1", "T2")
