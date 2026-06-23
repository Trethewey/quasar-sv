"""Tests for the artefact-mask + rescue inference pipeline."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.annotate import annotate_calls
from quasarsv.model import FusionCall
from quasarsv.qc import (
    flag_builtin_artefact_loci,
    flag_short_range_intrachr,
    ASHM_TARGETS,
)
from quasarsv.rescue import RescueConfig, rescue_ig_driver_pairs


def _fc(sample, ca, pa, cb, pb, gene_a="", gene_b="", tier="T2",
        sr=10, pe=0, callers=None):
    return FusionCall(
        sample=sample, fusion_id=f"{sample}__{ca}_{pa}_{cb}_{pb}",
        chrom_a=ca, pos_a=pa, strand_a="+",
        chrom_b=cb, pos_b=pb, strand_b="+",
        sv_type="BND", gene_a=gene_a, gene_b=gene_b,
        tier=tier, split_reads=sr, discordant_pairs=pe,
        callers_supporting=callers or ["forge_scan"], n_callers=1,
    )


def test_builtin_artefact_mask_flags_chr2_32916():
    # Inside artefact window
    c = _fc("S1", "3", 187_700_000, "2", 32_916_300, gene_a="BCL6")
    flag_builtin_artefact_loci([c])
    assert "builtin_artefact_locus" in c.qc_flags
    # Outside window
    c2 = _fc("S1", "3", 187_700_000, "2", 88_900_000, gene_a="BCL6", gene_b="IGK")
    flag_builtin_artefact_loci([c2])
    assert "builtin_artefact_locus" not in c2.qc_flags


def test_short_range_skips_ashm_targets():
    # Short-range intrachr at BCL6 — should NOT be flagged
    bcl6 = _fc("S1", "3", 187_700_000, "3", 187_700_100,
               gene_a="BCL6", gene_b="BCL6")
    flag_short_range_intrachr([bcl6])
    assert "short_range" not in bcl6.qc_flags
    # Same at a non-aSHM gene — should be flagged
    other = _fc("S1", "5", 100_000, "5", 100_080,
                gene_a="OTHER", gene_b="OTHER")
    flag_short_range_intrachr([other])
    assert "short_range" in other.qc_flags


def test_rescue_inference_emits_canonical_pair_with_high_support():
    art_chrom, art_pos = "2", 32_916_300
    # Driver BCL6 with strong artefact signal
    bcl6_to_art = _fc("S1", "3", 187_700_000, art_chrom, art_pos,
                      gene_a="BCL6", gene_b="", sr=400, pe=200)
    # IGH with strong artefact signal
    igh_to_art = _fc("S1", "14", 106_500_000, art_chrom, art_pos,
                     gene_a="IGH", gene_b="", sr=200, pe=100)
    bcl6_to_art.qc_flags.append("builtin_artefact_locus")
    igh_to_art.qc_flags.append("builtin_artefact_locus")
    calls = [bcl6_to_art, igh_to_art]
    rescue_ig_driver_pairs(calls, cfg=RescueConfig(min_artefact_sr_per_side=50))
    # At least one synthetic call should be emitted joining BCL6 + IGH
    rescued = [c for c in calls if "inferred_via_artefact_rescue" in c.qc_flags]
    assert any({c.gene_a, c.gene_b} == {"BCL6", "IGH"} for c in rescued)
    bcl6_igh = [c for c in rescued if {c.gene_a, c.gene_b} == {"BCL6", "IGH"}][0]
    # Canonical partner with known table → tier T1
    assert bcl6_igh.tier == "T1"
    assert bcl6_igh.known_partner is True
    assert "DLBCL" in bcl6_igh.known_partner_source


def test_rescue_caps_emissions_per_sample():
    art_pos = 32_916_300
    calls = []
    # Many drivers + many IGs → rescue could otherwise emit a fan-out
    for d, ch in [("BCL6", "3"), ("BCL2", "18"), ("CCND1", "11"),
                  ("MYC", "8"), ("MALT1", "18"), ("PAX5", "9")]:
        c = _fc("S1", ch, 50_000_000, "2", art_pos,
                gene_a=d, gene_b="", sr=300, pe=100)
        c.qc_flags.append("builtin_artefact_locus")
        calls.append(c)
    for ig, ch in [("IGH", "14"), ("IGK", "2"), ("IGL", "22"),
                   ("TRA", "14"), ("TRB", "7")]:
        c = _fc("S1", ch, 60_000_000, "2", art_pos,
                gene_a=ig, gene_b="", sr=300, pe=100)
        c.qc_flags.append("builtin_artefact_locus")
        calls.append(c)
    rescue_ig_driver_pairs(calls, cfg=RescueConfig(max_pairs_per_sample=4))
    rescued = [c for c in calls if "inferred_via_artefact_rescue" in c.qc_flags]
    assert len(rescued) <= 4


def test_rescue_bcell_lineage_excludes_tcell_receptors():
    """B-cell lineage prior must drop TRA/TRB even when their per-gene SR is high."""
    art_chrom, art_pos = "2", 32_916_300
    bcl6 = _fc("S1", "3", 187_700_000, art_chrom, art_pos, gene_a="BCL6", sr=400, pe=100)
    bcl6.qc_flags.append("builtin_artefact_locus")
    igh = _fc("S1", "14", 106_500_000, art_chrom, art_pos, gene_a="IGH", sr=100, pe=50)
    igh.qc_flags.append("builtin_artefact_locus")
    tra = _fc("S1", "14", 22_200_000, art_chrom, art_pos, gene_a="TRA", sr=2000, pe=500)
    tra.qc_flags.append("builtin_artefact_locus")
    calls = [bcl6, igh, tra]
    rescue_ig_driver_pairs(calls, cfg=RescueConfig(lineage="B",
                                                    min_artefact_sr_per_side=30))
    rescued = [c for c in calls if "inferred_via_artefact_rescue" in c.qc_flags]
    igs_emitted = {c.gene_b for c in rescued}
    assert "TRA" not in igs_emitted    # T-cell receptors excluded under B-cell prior
    assert "IGH" in igs_emitted        # IGH retained despite lower SR (canonical partner)


def test_rescue_emits_canonical_alternatives_with_ambiguous_flag():
    """When multiple canonical IGs have artefact signal, rescue emits all,
    flagging the lower-scoring ones with ``ig_partner_ambiguous``.
    Reproduces the PMBL Karpas-1106P pattern: BCL6 has BOTH IGL (high SR) and
    IGH (low SR) as artefact-mediated candidates; both are canonical t-band
    partners; both must be emitted so the truth (IGH) is visible."""
    art_chrom, art_pos = "2", 32_916_300
    bcl6 = _fc("S1", "3", 187_700_000, art_chrom, art_pos, gene_a="BCL6",
               sr=3000, pe=1000)
    bcl6.qc_flags.append("builtin_artefact_locus")
    # IGL with much HIGHER per-gene SR than IGH — mimics PMBL signal distribution
    igh = _fc("S1", "14", 106_500_000, art_chrom, art_pos, gene_a="IGH",
              sr=120, pe=80)
    igh.qc_flags.append("builtin_artefact_locus")
    igl = _fc("S1", "22", 22_300_000, art_chrom, art_pos, gene_a="IGL",
              sr=3300, pe=1200)
    igl.qc_flags.append("builtin_artefact_locus")
    igk = _fc("S1", "2", 89_500_000, art_chrom, art_pos, gene_a="IGK",
              sr=1400, pe=800)
    igk.qc_flags.append("builtin_artefact_locus")
    calls = [bcl6, igh, igl, igk]
    rescue_ig_driver_pairs(calls, cfg=RescueConfig(
        lineage="B",
        emit_canonical_alternatives=True,
        max_canonical_igs_per_driver=3,
        max_pairs_per_sample=12,
        min_artefact_sr_per_side=30,
    ))
    rescued = [c for c in calls if "inferred_via_artefact_rescue" in c.qc_flags]
    bcl6_pairs = {c.gene_b: c for c in rescued if c.gene_a == "BCL6"}
    # All three canonical IG partners must appear so the truth is recoverable
    assert "IGH" in bcl6_pairs, "IGH must be retained as a canonical alternative"
    assert "IGL" in bcl6_pairs
    assert "IGK" in bcl6_pairs
    # Top per-gene SR → primary (no ambiguity flag); the others tagged ambiguous
    assert "ig_partner_ambiguous" not in bcl6_pairs["IGL"].qc_flags
    assert "ig_partner_ambiguous" in bcl6_pairs["IGH"].qc_flags
    assert "ig_partner_ambiguous" in bcl6_pairs["IGK"].qc_flags
    # All three are canonical t-band partners → T1
    for ig in ("IGH", "IGL", "IGK"):
        assert bcl6_pairs[ig].tier == "T1"
        assert bcl6_pairs[ig].known_partner is True


def test_rescue_canonical_alternative_kept_when_other_driver_drops_it():
    """If the standard ratio_keep filter would drop IGH (because its per-gene
    SR is far below the top IG), the canonical-keep rule must rescue it."""
    art_chrom, art_pos = "2", 32_916_300
    bcl6 = _fc("S1", "3", 187_700_000, art_chrom, art_pos, gene_a="BCL6",
               sr=3000, pe=1000)
    bcl6.qc_flags.append("builtin_artefact_locus")
    # IGH SR is 3 % of IGL — below ratio_keep=0.20 cutoff
    igh = _fc("S1", "14", 106_500_000, art_chrom, art_pos, gene_a="IGH",
              sr=120, pe=80)
    igh.qc_flags.append("builtin_artefact_locus")
    igl = _fc("S1", "22", 22_300_000, art_chrom, art_pos, gene_a="IGL",
              sr=4000, pe=1500)
    igl.qc_flags.append("builtin_artefact_locus")
    calls = [bcl6, igh, igl]
    rescue_ig_driver_pairs(calls, cfg=RescueConfig(
        lineage="B", emit_canonical_alternatives=True,
        ratio_keep=0.20, min_artefact_sr_per_side=30,
    ))
    rescued = [c for c in calls if "inferred_via_artefact_rescue" in c.qc_flags]
    bcl6_igh = [c for c in rescued if c.gene_a == "BCL6" and c.gene_b == "IGH"]
    assert bcl6_igh, "BCL6-IGH must be emitted via canonical-keep, not dropped by ratio_keep"


def test_rescue_tcell_lineage_excludes_bcell_igs():
    """Symmetry: T-cell lineage prior keeps TRA/TRB and drops IGH/IGK/IGL."""
    art_chrom, art_pos = "2", 32_916_300
    alk = _fc("S2", "2", 29_500_000, art_chrom, art_pos, gene_a="ALK", sr=400, pe=200)
    alk.qc_flags.append("builtin_artefact_locus")
    tra = _fc("S2", "14", 22_200_000, art_chrom, art_pos, gene_a="TRA", sr=500, pe=300)
    tra.qc_flags.append("builtin_artefact_locus")
    igh = _fc("S2", "14", 106_500_000, art_chrom, art_pos, gene_a="IGH", sr=2000, pe=1000)
    igh.qc_flags.append("builtin_artefact_locus")
    calls = [alk, tra, igh]
    rescue_ig_driver_pairs(
        calls,
        cfg=RescueConfig(lineage="T", min_artefact_sr_per_side=30),
    )
    rescued = [c for c in calls if "inferred_via_artefact_rescue" in c.qc_flags]
    igs_emitted = {c.gene_b for c in rescued}
    assert "IGH" not in igs_emitted    # B-cell IG dropped under T-cell prior
    # ALK has no canonical TR partner in known_partners.tsv, so the only pair is
    # non-canonical ALK-TRA; lineage still routes it through the T-cell IG set.
    assert "TRA" in igs_emitted


def test_rescue_per_sample_lineage_override():
    """sample_lineage map must override the default lineage per sample."""
    art_chrom, art_pos = "2", 32_916_300
    # B-cell sample S1 with TRA signal — TRA should be excluded
    bcl6_s1 = _fc("S1", "3", 187_700_000, art_chrom, art_pos, gene_a="BCL6", sr=500, pe=100)
    bcl6_s1.qc_flags.append("builtin_artefact_locus")
    tra_s1 = _fc("S1", "14", 22_200_000, art_chrom, art_pos, gene_a="TRA", sr=400, pe=100)
    tra_s1.qc_flags.append("builtin_artefact_locus")
    igh_s1 = _fc("S1", "14", 106_500_000, art_chrom, art_pos, gene_a="IGH", sr=200, pe=100)
    igh_s1.qc_flags.append("builtin_artefact_locus")
    # T-cell sample S2 with IGH signal — IGH should be excluded
    alk_s2 = _fc("S2", "2", 29_500_000, art_chrom, art_pos, gene_a="ALK", sr=500, pe=100)
    alk_s2.qc_flags.append("builtin_artefact_locus")
    tra_s2 = _fc("S2", "14", 22_200_000, art_chrom, art_pos, gene_a="TRA", sr=400, pe=100)
    tra_s2.qc_flags.append("builtin_artefact_locus")
    igh_s2 = _fc("S2", "14", 106_500_000, art_chrom, art_pos, gene_a="IGH", sr=2000, pe=1000)
    igh_s2.qc_flags.append("builtin_artefact_locus")
    calls = [bcl6_s1, tra_s1, igh_s1, alk_s2, tra_s2, igh_s2]
    rescue_ig_driver_pairs(
        calls,
        cfg=RescueConfig(min_artefact_sr_per_side=30),
        sample_lineage={"S1": "B", "S2": "T"},
    )
    s1_rescued = [c for c in calls if c.sample == "S1" and "inferred_via_artefact_rescue" in c.qc_flags]
    s2_rescued = [c for c in calls if c.sample == "S2" and "inferred_via_artefact_rescue" in c.qc_flags]
    s1_igs = {c.gene_b for c in s1_rescued}
    s2_igs = {c.gene_b for c in s2_rescued}
    assert "TRA" not in s1_igs and "IGH" in s1_igs
    assert "IGH" not in s2_igs and "TRA" in s2_igs


def test_rescue_suppresses_canonical_fanout_from_weak_drivers():
    """When the top canonical driver has strong artefact signal, weak
    secondary drivers that just happen to be canonical partners of the same
    shared IG must NOT be emitted as additional T1 pairs. (Otherwise BCL2-IGH,
    CCND1-IGH, MALT1-IGH all fan out from one shared low IGH signal.)"""
    art_chrom, art_pos = "2", 32_916_300
    # Top canonical driver: BCL6 with strong signal
    bcl6 = _fc("S1", "3", 187_700_000, art_chrom, art_pos, gene_a="BCL6", sr=3000, pe=1200)
    bcl6.qc_flags.append("builtin_artefact_locus")
    # Weak passenger drivers also at the artefact (10 % of BCL6 signal)
    bcl2 = _fc("S1", "18", 63_200_000, art_chrom, art_pos, gene_a="BCL2", sr=300, pe=80)
    bcl2.qc_flags.append("builtin_artefact_locus")
    ccnd1 = _fc("S1", "11", 69_645_000, art_chrom, art_pos, gene_a="CCND1", sr=300, pe=80)
    ccnd1.qc_flags.append("builtin_artefact_locus")
    malt1 = _fc("S1", "18", 58_700_000, art_chrom, art_pos, gene_a="MALT1", sr=300, pe=80)
    malt1.qc_flags.append("builtin_artefact_locus")
    # Shared IGH signal — low
    igh = _fc("S1", "14", 106_500_000, art_chrom, art_pos, gene_a="IGH", sr=100, pe=60)
    igh.qc_flags.append("builtin_artefact_locus")
    # Strong shared IGL — IGL is BCL6's dominant artefact-mediated partner
    igl = _fc("S1", "22", 22_300_000, art_chrom, art_pos, gene_a="IGL", sr=3300, pe=1200)
    igl.qc_flags.append("builtin_artefact_locus")
    calls = [bcl6, bcl2, ccnd1, malt1, igh, igl]
    rescue_ig_driver_pairs(calls, cfg=RescueConfig(
        lineage="B", emit_canonical_alternatives=True,
        noncanonical_fanout_ratio=0.20,
        min_artefact_sr_per_side=30,
    ))
    rescued = [c for c in calls if "inferred_via_artefact_rescue" in c.qc_flags]
    pairs = {(c.gene_a, c.gene_b) for c in rescued}
    # BCL6 (top canonical) emits IGL/IGH primary + alts
    assert ("BCL6", "IGL") in pairs
    assert ("BCL6", "IGH") in pairs
    # BCL2-IGH, CCND1-IGH, MALT1-IGH share the same weak IGH signal — all
    # have score 100 << 0.20 × 3000 = 600 → suppressed
    assert ("BCL2", "IGH") not in pairs
    assert ("CCND1", "IGH") not in pairs
    assert ("MALT1", "IGH") not in pairs


def test_metadata_cohort_to_lineage():
    from quasarsv.metadata import cohort_to_lineage
    assert cohort_to_lineage("PMBL") == "B"
    assert cohort_to_lineage("DLBCL") == "B"
    assert cohort_to_lineage("FL") == "B"
    assert cohort_to_lineage("ATLL") == "T"
    assert cohort_to_lineage("PTCL") == "T"
    assert cohort_to_lineage("T-cell NHL") == "T"
    assert cohort_to_lineage("") == "any"
    assert cohort_to_lineage("Unknown") == "any"
