"""Tests for library_stats — focused on the MAPQ weighting + adaptive-insert math.

The pysam-driven CRAM read step is excluded from CI tests (no fixture BAM
shipped). Those paths are exercised indirectly by the WGS smoke tests.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.scanners.library_stats import (
    LibraryStats, mapq_weight, _percentile,
)


def test_mapq_weight_thresholds():
    # MAPQ 0 → 0.0
    assert mapq_weight(0) == 0.0
    # MAPQ 30 → 1.0 (default full_weight_mapq)
    assert mapq_weight(30) == 1.0
    # MAPQ 60 → 1.0 (saturates)
    assert mapq_weight(60) == 1.0
    # MAPQ 15 → 0.5 (linear taper from 0 to 30)
    assert mapq_weight(15) == 0.5
    # Custom full-weight floor
    assert mapq_weight(10, full_weight_mapq=10) == 1.0
    assert mapq_weight(5, full_weight_mapq=10) == 0.5


def test_adaptive_insert_threshold_uses_median_plus_5mad():
    lib = LibraryStats(insert_median=350.0, insert_mad=40.0)
    # median + 5 * MAD = 350 + 200 = 550
    assert lib.discordant_min_distance == 550


def test_adaptive_insert_falls_back_when_no_inserts():
    """No insert sizes observed → default 10 kb."""
    lib = LibraryStats(insert_median=0.0, insert_mad=0.0)
    assert lib.discordant_min_distance == 10_000


def test_library_stats_round_trips_json():
    orig = LibraryStats(
        n_reads_sampled=12345,
        insert_median=380.0, insert_mad=42.0,
        insert_p95=580.0, insert_p99=720.0,
        softclip_median=0.0, softclip_p95=8.0,
        mapq_median=60.0, mapq_p10=20.0,
        pct_supplementary=0.01, pct_duplicate=0.05,
    )
    s = orig.to_json()
    back = LibraryStats.from_json(s)
    assert back.insert_median == 380.0
    assert back.discordant_min_distance == orig.discordant_min_distance
    assert back.pct_duplicate == 0.05


def test_library_stats_save_load(tmp_path):
    p = tmp_path / "stats.json"
    lib = LibraryStats(insert_median=400.0, insert_mad=50.0)
    lib.save(p)
    assert p.exists()
    back = LibraryStats.load(p)
    assert back.insert_median == 400.0


def test_percentile_helper():
    vals = sorted([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    assert _percentile(vals, 0.0) == 1
    assert _percentile(vals, 1.0) == 10
    # 0.5 of 0-indexed 0..9 = index 4.5 → 5.5
    assert _percentile(vals, 0.5) == 5.5
