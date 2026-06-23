"""Tests for the DBSCAN-style clustering in merge.py.

The DBSCAN clusterer replaces fixed-grid buckets with density-based union-find.
Key properties to verify:

1. Two close calls always cluster (the legacy bucket approach could split them
   across bucket boundaries — the DBSCAN approach never does).
2. Three calls forming a chain within eps each-pairwise-step transitively cluster.
3. Calls with different chroms / strands never cluster, even if positions match.
4. ``min_samples=2`` drops singletons.
5. Behaviour matches legacy bucket clustering on isolated points.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.merge import MergeConfig, cluster_calls, _cluster_bucket, _cluster_dbscan
from quasarsv.model import BreakpointCall, Evidence


def _bp(chrom_a="14", pos_a=106_000_000, strand_a="+",
        chrom_b="18", pos_b=63_200_000, strand_b="-",
        caller="manta", sr=5, pe=10):
    return BreakpointCall(
        sample="S1",
        chrom_a=chrom_a, pos_a=pos_a, strand_a=strand_a,
        chrom_b=chrom_b, pos_b=pos_b, strand_b=strand_b,
        sv_type="BND",
        evidence=Evidence(caller=caller, split_reads=sr, discordant_pairs=pe,
                          filter_pass=True, precise=True, raw_qual=20.0),
        record_id=f"r_{pos_a}_{pos_b}_{caller}",
    )


def test_dbscan_close_calls_cluster():
    """Two calls within eps in both dimensions form one cluster."""
    calls = [_bp(pos_a=106_000_000, pos_b=63_200_000, caller="manta"),
             _bp(pos_a=106_000_100, pos_b=63_200_100, caller="delly")]
    cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=250))
    assert len(cands) == 1
    assert set(cands[0].member_callers) == {"manta", "delly"}


def test_dbscan_eliminates_bucket_boundary_artefact():
    """Two calls straddling a bucket boundary still cluster in DBSCAN —
    where the legacy bucket approach can drop the link."""
    # Bucket boundary at pos // 250 — calls at 249 and 251 fall in different
    # buckets but are only 2 bp apart.
    calls = [_bp(pos_a=106_000_249, pos_b=63_200_249, caller="manta"),
             _bp(pos_a=106_000_251, pos_b=63_200_251, caller="delly")]
    cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=250))
    # DBSCAN should join: they're 2 bp apart in both dims
    assert len(cands) == 1
    assert set(cands[0].member_callers) == {"manta", "delly"}


def test_dbscan_transitive_chain():
    """A → B → C chain where A-C exceeds eps but A-B and B-C are within eps
    transitively forms one cluster."""
    calls = [_bp(pos_a=106_000_000, pos_b=63_200_000, caller="m1"),
             _bp(pos_a=106_000_200, pos_b=63_200_200, caller="m2"),
             _bp(pos_a=106_000_400, pos_b=63_200_400, caller="m3")]
    # eps=250: m1-m2 dist=200 (within), m2-m3 dist=200 (within), m1-m3 dist=400 (NOT)
    cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=250))
    assert len(cands) == 1     # all three transitively in same cluster
    assert set(cands[0].member_callers) == {"m1", "m2", "m3"}


def test_dbscan_different_chroms_dont_cluster():
    calls = [_bp(chrom_a="14", pos_a=100, chrom_b="18", pos_b=100, caller="manta"),
             _bp(chrom_a="14", pos_a=100, chrom_b="19", pos_b=100, caller="delly")]
    cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=10_000))
    assert len(cands) == 2


def test_dbscan_different_strands_dont_cluster_when_required():
    calls = [_bp(strand_a="+", strand_b="-", caller="manta"),
             _bp(strand_a="-", strand_b="+", caller="delly")]
    cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=10_000,
                                                same_strand_required=True))
    assert len(cands) == 2


def test_dbscan_different_strands_cluster_when_not_required():
    calls = [_bp(strand_a="+", strand_b="-", caller="manta"),
             _bp(strand_a="-", strand_b="+", caller="delly")]
    cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=10_000,
                                                same_strand_required=False))
    assert len(cands) == 1


def test_dbscan_min_samples_drops_singletons():
    calls = [_bp(pos_a=1_000_000, pos_b=2_000_000, caller="manta"),
             _bp(pos_a=8_000_000, pos_b=9_000_000, caller="delly")]
    cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=250,
                                                dbscan_min_samples=2))
    assert cands == []     # two isolated singletons, both dropped


def test_dbscan_default_keeps_singletons():
    calls = [_bp(pos_a=1_000_000, pos_b=2_000_000, caller="manta")]
    cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=250,
                                                dbscan_min_samples=1))
    assert len(cands) == 1


def test_dbscan_handles_empty_input():
    assert _cluster_dbscan([], MergeConfig()) == []


def test_top_level_cluster_calls_uses_dbscan_by_default():
    """The public API picks DBSCAN by default per docs/precision_techniques.md #7."""
    cfg = MergeConfig()
    assert cfg.clustering == "dbscan"
    calls = [_bp(pos_a=106_000_249, pos_b=63_200_249, caller="manta"),
             _bp(pos_a=106_000_251, pos_b=63_200_251, caller="delly")]
    cands = cluster_calls(calls, cfg)
    # DBSCAN joins; bucket clustering would split — verify we get one cluster.
    assert len(cands) == 1


def test_cluster_calls_bucket_mode_still_works():
    """Legacy bucket mode is still accessible for regression / debug."""
    cfg = MergeConfig(clustering="bucket")
    calls = [_bp(pos_a=106_000_000, pos_b=63_200_000, caller="manta"),
             _bp(pos_a=106_000_100, pos_b=63_200_100, caller="delly")]
    cands = cluster_calls(calls, cfg)
    assert len(cands) == 1


def test_dbscan_matches_bucket_on_isolated_clusters():
    """When all clusters sit well inside bucket boundaries, DBSCAN and
    bucket should produce the same partitioning."""
    calls = [
        # Cluster 1
        _bp(pos_a=10_000, pos_b=20_000, caller="manta"),
        _bp(pos_a=10_050, pos_b=20_050, caller="delly"),
        # Cluster 2 (far away)
        _bp(pos_a=80_000, pos_b=90_000, caller="manta"),
        _bp(pos_a=80_050, pos_b=90_050, caller="delly"),
    ]
    db_cands = _cluster_dbscan(calls, MergeConfig(pos_tolerance=250))
    bu_cands = _cluster_bucket(calls, MergeConfig(pos_tolerance=250))
    assert len(db_cands) == 2
    assert len(bu_cands) == 2
    db_partition = sorted(sorted(c.member_callers) for c in db_cands)
    bu_partition = sorted(sorted(c.member_callers) for c in bu_cands)
    assert db_partition == bu_partition
