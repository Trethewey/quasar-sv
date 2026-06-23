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


def test_truthset_ambiguous_alts_default_not_fp():
    """Rescue's intentional ambiguous candidates must NOT degrade precision."""
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [
        _fc("S1", "BCL6", "IGL", tier="T1"),                                      # primary mis-pick
        _fc("S1", "BCL6", "IGK", tier="T1", flags=["ig_partner_ambiguous"]),      # ambiguous alt
        _fc("S1", "BCL6", "IGH", tier="T1", flags=["ig_partner_ambiguous"]),      # ambiguous alt — TRUTH
    ]
    bench = score_calls_against_truth(calls, truth)
    assert bench.tp == 1       # truth IS in the call set (the ambiguous alt)
    # Only BCL6-IGL counts as FP (BCL6-IGK is ambiguous-flagged → exempt by default)
    assert bench.fp == 1


def test_truthset_ambiguous_alts_strict_mode():
    """Strict mode counts every non-truth T1 call as FP, ambiguous or not."""
    truth = [TruthEntry(sample_id="S1", gene_a="BCL6", gene_b="IGH",
                        truth_class="confirmed")]
    calls = [
        _fc("S1", "BCL6", "IGL", tier="T1"),
        _fc("S1", "BCL6", "IGK", tier="T1", flags=["ig_partner_ambiguous"]),
        _fc("S1", "BCL6", "IGH", tier="T1", flags=["ig_partner_ambiguous"]),
    ]
    bench = score_calls_against_truth(calls, truth, ambiguous_alts_count_as_fp=True)
    assert bench.tp == 1
    assert bench.fp == 2       # BCL6-IGL + BCL6-IGK (BCL6-IGH matches truth)


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
    # Karpas-1106P PMBL truth should be present
    karpas = [t for t in truth if "Karpas1106P" in t.sample_id]
    assert karpas
    assert karpas[0].pair == frozenset({"BCL6", "IGH"})
    assert karpas[0].is_positive
