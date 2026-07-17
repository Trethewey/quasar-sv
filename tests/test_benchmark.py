"""Tests for the truth-set scorer."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.benchmark import (
    TruthEntry, score_calls_against_truth, load_truth_set, builtin_truth_path,
)
from quasarsv.model import FusionCall


def _fc(sample, ga, gb, tier="T1", flags=None):
    return FusionCall(
        sample=sample, fusion_id=f"{sample}__{ga}_{gb}",
        chrom_a="14", pos_a=106_000_000, strand_a="+",
        chrom_b="3", pos_b=187_500_000, strand_b="+",
        sv_type="BND", gene_a=ga, gene_b=gb,
        tier=tier, qc_flags=list(flags or []),
    )


def test_truthset_matches_canonical_pair_order_insensitive():
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [_fc("S1", "IGH", "BCL6", tier="T1")]    # reversed order
    bench = score_calls_against_truth(calls, truth)
    assert bench.tp == 1
    assert bench.fn == 0
    assert bench.recall == 1.0


def test_truthset_negative_control_t1_calls_are_fp():
    truth = [TruthEntry(sample_id="S1", truth_class="none_expected")]
    calls = [_fc("S1", "BCL6", "IGH", tier="T1")]
    bench = score_calls_against_truth(calls, truth)
    assert bench.tp == 0
    assert bench.fn == 0
    assert bench.fp == 1
    assert bench.per_sample["S1"].is_negative_control


def test_truthset_ambiguous_alts_count_as_fp_by_default():
    """Hedged alternates must cost precision.

    Emitting IGH, IGK *and* IGL for one driver and being scored only on the one
    that happens to hit the truth is not precision — it is a spread bet. The
    default must charge for the misses.
    """
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [
        _fc("S1", "BCL6", "IGL", tier="T1"),                                      # primary mis-pick
        _fc("S1", "BCL6", "IGK", tier="T1", flags=["ig_partner_ambiguous"]),      # hedged alt
        _fc("S1", "BCL6", "IGH", tier="T1", flags=["ig_partner_ambiguous"]),      # hedged alt — TRUTH
    ]
    bench = score_calls_against_truth(calls, truth)
    assert bench.tp == 1
    assert bench.fp == 2       # BCL6-IGL + BCL6-IGK both charged


def test_truthset_ambiguous_alts_lenient_opt_out():
    """The old lenient behaviour survives only as an explicit opt-out."""
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [
        _fc("S1", "BCL6", "IGL", tier="T1"),
        _fc("S1", "BCL6", "IGK", tier="T1", flags=["ig_partner_ambiguous"]),
        _fc("S1", "BCL6", "IGH", tier="T1", flags=["ig_partner_ambiguous"]),
    ]
    bench = score_calls_against_truth(calls, truth, ambiguous_alts_count_as_fp=False)
    assert bench.tp == 1
    assert bench.fp == 1       # only the unflagged BCL6-IGL


def test_driver_only_does_not_match_driver_ig_truth_by_default():
    """"Found the driver, couldn't resolve the partner" is not a detected fusion."""
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [_fc("S1", "BCL6", "", tier="T1")]
    assert score_calls_against_truth(calls, truth).tp == 0
    # Only a deliberately lenient side-metric may count it.
    lenient = score_calls_against_truth(
        calls, truth, relax_canonical_ig_partner=True, match_driver_only=True)
    assert lenient.tp == 1


def test_detection_vs_lookup_split():
    """A TP with no read-level junction is a lookup, and must be reported as one."""
    truth = [
        TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH", truth_class="confirmed"),
        TruthEntry(sample_id="S2", gene_a="BCL2", gene_b="IGH", truth_class="confirmed"),
    ]
    calls = [
        _fc("S1", "BCL6", "IGH", tier="T1"),      # right answer, NO junction in reads
        _fc("S2", "BCL2", "IGH", tier="T1"),      # right answer, junction present
    ]
    support = {("S2", frozenset({"BCL2", "IGH"}))}
    bench = score_calls_against_truth(calls, truth, junction_support=support)
    assert bench.tp == 2
    assert bench.tp_detected == 1
    assert bench.tp_lookup_only == 1


def test_truthset_missing_canonical_is_fn():
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [_fc("S1", "ALK", "IGL", tier="T2")]     # truth not in calls
    bench = score_calls_against_truth(calls, truth)
    assert bench.tp == 0
    assert bench.fn == 1
    assert bench.recall == 0.0


def test_truthset_double_hit_dlbcl():
    """A sample with two canonical translocations counts both."""
    truth = [
        TruthEntry(sample_id="DH1", gene_a="IGH", gene_b="BCL2", truth_class="confirmed"),
        TruthEntry(sample_id="DH1", gene_a="MYC", gene_b="IGH", truth_class="confirmed"),
    ]
    calls = [
        _fc("DH1", "IGH", "BCL2", tier="T1"),
        _fc("DH1", "MYC", "IGH", tier="T1"),
    ]
    bench = score_calls_against_truth(calls, truth)
    assert bench.tp == 2
    assert bench.fn == 0


def test_truthset_match_tier_reports_best():
    """When the same truth pair appears at multiple tiers, report the best."""
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [
        _fc("S1", "BCL6", "IGH", tier="T2"),
        _fc("S1", "BCL6", "IGH", tier="T1"),
    ]
    bench = score_calls_against_truth(calls, truth)
    assert bench.per_sample["S1"].match_tier[frozenset({"BCL6", "IGH"})] == "T1"


def test_builtin_truth_loads():
    """The packaged cohort_truth.tsv must be parseable."""
    truth = load_truth_set(builtin_truth_path())
    assert len(truth) > 0
    karpas = [t for t in truth if "Karpas1106P" in t.sample_id]
    assert karpas


def test_karpas_and_u2940_are_negative_controls():
    """Karpas-1106P and U2940 must NOT carry a BCL6-IGH truth.

    Break-apart FISH over BCL6 and IG-H/K/L is germline in both lines (Dai 2015,
    PMID:26599546), and the read-level oracle finds zero BCL6-IGH reads of either
    kind. The t(3;14) previously recorded for them is refuted; they are negative
    controls. This test exists because the old truth entry, paired with a
    canonical-partner prior, was what let a fabricated BCL6-IGH score as a true
    positive.
    """
    truth = load_truth_set(builtin_truth_path())
    for sid in ("ERR9188549_Karpas1106P", "ERR9128954_U2940"):
        rows = [t for t in truth if t.sample_id == sid]
        assert rows, f"{sid} missing from truth set"
        assert all(t.is_negative for t in rows), f"{sid} must be a negative control"
        assert not any(t.pair == frozenset({"BCL6", "IGH"}) and t.is_positive
                       for t in rows)


def test_md903_keeps_its_real_bcl6_igh():
    """The genuine BCL6-IGH must survive the correction.

    MD903 is the cohort's positive control for t(3;14): its partner is stated
    explicitly in the primary literature and the read-level oracle sees it at
    PE=16. It proves a real BCL6-IGH is detectable, so the correction above
    removes fabrications without removing the event class.
    """
    truth = load_truth_set(builtin_truth_path())
    rows = [t for t in truth if t.sample_id == "SRR1236472_MD903_DLBCL_cell_line"]
    assert any(t.pair == frozenset({"BCL6", "IGH"}) and t.is_positive for t in rows)


def test_disputed_truth_is_quarantined_not_scored():
    """A disputed row is neither a positive, a negative, nor an FP magnet."""
    truth = [
        TruthEntry(sample_id="S1", gene_a="MYC", gene_b="IGL", truth_class="disputed"),
    ]
    calls = [_fc("S1", "MYC", "IGL", tier="T1"), _fc("S1", "ALK", "IGH", tier="T1")]
    bench = score_calls_against_truth(calls, truth)
    # The sample is unscoreable: it must not appear at all, and in particular
    # must not turn every call into a false positive via an empty expected set.
    assert bench.tp == 0 and bench.fn == 0 and bench.fp == 0
    assert "S1" not in bench.per_sample


def test_fp_charged_at_same_tiers_as_tp_by_default():
    """Free recall guard: a T2 call that hits truth is a TP, so a T2 call that
    misses must be an FP. Otherwise volume at T2 is a pure win."""
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [
        _fc("S1", "BCL6", "IGH", tier="T2"),     # TP at T2
        _fc("S1", "ALK", "IGK", tier="T2"),      # miss at T2 -> must cost
    ]
    bench = score_calls_against_truth(calls, truth)
    assert bench.tp == 1
    assert bench.fp == 1
    # The clinically-actionable view charges (and credits) T1 only.
    t1_only = score_calls_against_truth(calls, truth, fp_tiers=("T1",))
    assert t1_only.fp == 0
