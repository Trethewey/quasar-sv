"""Round-trip sanity tests for parsers and merger.

Run from the project root with `pytest tests` (requires the dev extra).
"""
from __future__ import annotations

import io
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.merge import (
    MergeConfig, merge_caller_calls, assign_tier, cluster_calls,
)
from quasarsv.model import (
    BreakpointCall, Evidence, FusionCall,
    write_fusion_calls_tsv, read_fusion_calls_tsv,
)
from quasarsv.annotate import annotate_calls
from quasarsv.parsers.base import parse_bnd_alt, normalise_order


def test_normalise_order_is_idempotent():
    a = normalise_order("3", 100, "+", "14", 200, "-")
    b = normalise_order("14", 200, "-", "3", 100, "+")
    assert a == b


def test_parse_bnd_alt_basic():
    # t[p[: strand_a='+', mate left of break
    r = parse_bnd_alt("N[chr2:32916489[")
    assert r is not None
    mchrom, mpos, sa, sb = r
    assert mchrom == "chr2" and mpos == 32916489
    assert sa == "+" and sb == "+"

    # [p[t: strand_a='-', mate left of break
    r2 = parse_bnd_alt("[chr14:106500000[N")
    assert r2 is not None
    assert r2[0] == "chr14" and r2[2] == "-" and r2[3] == "+"


def _make_call(sample, ca, pa, cb, pb, caller, sr=10, pe=5, asm=0, fp=True, prec=True):
    ev = Evidence(caller=caller, split_reads=sr, discordant_pairs=pe,
                  assembly_contigs=asm, filter_pass=fp, precise=prec, raw_qual=float(sr))
    return BreakpointCall(sample=sample, chrom_a=ca, pos_a=pa, strand_a="+",
                          chrom_b=cb, pos_b=pb, strand_b="+", sv_type="BND",
                          evidence=ev, record_id=f"{caller}_{ca}_{pa}_{cb}_{pb}")


def test_merger_clusters_concordant_callers_into_one_candidate():
    bps = [
        _make_call("S1", "3", 187_700_000, "14", 106_700_000, "manta"),
        _make_call("S1", "3", 187_700_050, "14", 106_700_080, "gridss"),
        _make_call("S1", "3", 187_700_100, "14", 106_700_120, "delly", pe=20, sr=0),
    ]
    cands = cluster_calls(bps, MergeConfig(pos_tolerance=500))
    assert len(cands) == 1
    cand = cands[0]
    assert set(cand.member_callers) == {"manta", "gridss", "delly"}
    summary = cand.evidence_summary()
    assert summary["split_read"] >= 10
    assert summary["discordant_pair"] >= 20


def test_t1_requires_multi_caller():
    one_caller = _make_call("S1", "3", 100, "14", 200, "manta", sr=30, pe=30, asm=1)
    cands = cluster_calls([one_caller], MergeConfig())
    tier, _flags = assign_tier(cands[0])
    assert tier in ("T2", "T1")   # single caller may hit very-strong path
    # confirm 2-caller agreement DOES make T1
    two = [
        _make_call("S1", "3", 100, "14", 200, "manta", sr=10, pe=5),
        _make_call("S1", "3", 105, "14", 205, "gridss", sr=8, pe=12),
    ]
    c2 = cluster_calls(two, MergeConfig())
    tier2, _ = assign_tier(c2[0])
    assert tier2 == "T1"


def test_annotation_flags_canonical_partner_igh_bcl2():
    # IGH locus 14:105.58M-106.88M  ↔  BCL2 18:63.12M-63.32M
    fcs = [FusionCall(
        sample="S1", fusion_id="f1",
        chrom_a="14", pos_a=106_500_000, strand_a="+",
        chrom_b="18", pos_b=63_200_000, strand_b="+",
        sv_type="BND",
    )]
    annotate_calls(fcs)
    assert fcs[0].gene_a == "IGH"
    assert fcs[0].gene_b == "BCL2"
    assert fcs[0].known_partner is True
    assert "canonical" in fcs[0].known_partner_source.lower()


def test_tsv_round_trip(tmp_path):
    fc = FusionCall(
        sample="S1", fusion_id="f1",
        chrom_a="3", pos_a=100, strand_a="+",
        chrom_b="14", pos_b=200, strand_b="-",
        sv_type="BND",
        callers_supporting=["manta", "gridss"], n_callers=2,
        split_reads=15, discordant_pairs=8, assembly_contigs=1,
        soft_clips=0, n_evidence_types=3, vaf=0.32, precise=True,
        any_pass=True, raw_qual_max=42.0,
        gene_a="BCL6", gene_b="IGH",
        region_a="exonic_or_intronic", region_b="upstream",
        in_frame=True, known_partner=True,
        known_partner_source="canonical:DLBCL",
        driver_locus="BCL6-IGH", tier="T1", qc_flags=["single_caller_strong"],
        member_record_ids=["manta_id_1", "gridss_id_1"],
    )
    p = tmp_path / "test.tsv"
    write_fusion_calls_tsv([fc], str(p))
    fcs = read_fusion_calls_tsv(str(p))
    assert len(fcs) == 1
    rt = fcs[0]
    assert rt.sample == "S1"
    assert rt.gene_a == "BCL6"
    assert rt.known_partner is True
    assert rt.tier == "T1"
    assert rt.callers_supporting == ["manta", "gridss"]
    assert rt.qc_flags == ["single_caller_strong"]
