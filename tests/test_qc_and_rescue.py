"""Tests for the artefact-mask pipeline and the refusal to fabricate partners.

The artefact channel at chr2:32,916 is adapter read-through and 2-colour poly-G
tails, emitted at a uniform rate by every locus in the genome (measured: ~220
split reads per 10k reads for rearranged drivers and unrearranged controls
alike). It therefore carries no information about which loci are partners.

These tests pin the resulting contract: artefact co-occurrence must never
produce a named fusion partner, at any tier, under any configuration.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.annotate import annotate_calls
from quasarsv.model import FusionCall
from quasarsv.qc import (
    flag_builtin_artefact_loci,
    flag_short_range_intrachr,
    demote_interchromosomal_without_discordant,
    ASHM_TARGETS,
)
from quasarsv.rescue import (
    RescueConfig,
    flag_artefact_masked_breakends,
    rescue_ig_driver_pairs,
)
from quasarsv.scanners.cram_scanner import is_noise_clip


def _fc(sample, ca, pa, cb, pb, gene_a="", gene_b="", tier="T2",
        sr=10, pe=0, callers=None):
    return FusionCall(
        sample=sample, fusion_id=f"{sample}__{ca}_{pa}_{cb}_{pb}",
        chrom_a=ca, pos_a=pa, strand_a="+",
        chrom_b=cb, pos_b=pb, strand_b="+",
        sv_type="BND", gene_a=gene_a, gene_b=gene_b,
        tier=tier, split_reads=sr, discordant_pairs=pe,
        callers_supporting=callers or ["quasar"], n_callers=1,
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


# ---- the fabrication must stay dead ----

def test_partner_inference_from_artefact_is_removed():
    """The old entry point must not silently resurrect."""
    with pytest.raises(NotImplementedError):
        rescue_ig_driver_pairs([], cfg=RescueConfig())


def test_removed_flags_are_rejected_not_silently_ignored():
    """A flag the pipeline no longer reads must ERROR, not be quietly accepted.

    --lineage / --metadata fed only the partner inference, and
    --chrom-sa-inference drove the coordinate-fabricating module. Leaving any of
    them parseable would let an operator believe they had set something. Hiding
    one behind argparse.SUPPRESS is not removal — it still parses.
    """
    from quasarsv.cli import main
    for flag, val in (("--lineage", "B"), ("--metadata", "m.xlsx"),
                      ("--chrom-sa-inference", None)):
        argv = ["call", "--sample", "S", "--bam", "x.cram",
                "--reference", "r.fa", "--output-dir", "o", flag]
        if val is not None:
            argv.append(val)
        # argparse exits 2 on an unrecognised argument; anything else means the
        # flag was accepted and silently ignored.
        with pytest.raises(SystemExit) as e:
            main(argv)
        assert e.value.code == 2, f"{flag} was accepted rather than rejected"


def test_lineage_prior_is_gone_not_merely_ignored():
    """The B/T prior must not survive as a control that does nothing.

    Its only consumer was the rescue's choice of which IG locus to NAME as a
    partner. With the naming removed, a retained `lineage` field (or a CLI
    --lineage flag) would advertise a prior the pipeline never reads — the same
    class of dishonesty as the fabrication itself.
    """
    assert not hasattr(RescueConfig(), "lineage")
    import inspect
    from quasarsv.qc import apply_default_qc
    params = inspect.signature(apply_default_qc).parameters
    assert "sample_lineage" not in params
    assert "lineage_default" not in params


def test_artefact_cooccurrence_never_invents_a_partner():
    """The exact Karpas-1106P pattern that used to yield a T1 BCL6-IGH call.

    A driver and several IG loci all dump reads into the artefact. No read joins
    BCL6 to any of them. Nothing may be synthesised — not IGH (the canonical
    prior), not IGL (the highest count), not anything.
    """
    art_chrom, art_pos = "2", 32_916_300
    bcl6 = _fc("S1", "3", 187_700_000, art_chrom, art_pos, gene_a="BCL6",
               sr=3000, pe=1000)
    igh = _fc("S1", "14", 106_500_000, art_chrom, art_pos, gene_a="IGH",
              sr=120, pe=80)
    igl = _fc("S1", "22", 22_300_000, art_chrom, art_pos, gene_a="IGL",
              sr=3300, pe=1200)
    igk = _fc("S1", "2", 89_500_000, art_chrom, art_pos, gene_a="IGK",
              sr=1400, pe=800)
    calls = [bcl6, igh, igl, igk]
    for c in calls:
        c.qc_flags.append("builtin_artefact_locus")

    n_before = len(calls)
    flag_artefact_masked_breakends(calls, cfg=RescueConfig())

    assert len(calls) == n_before, "no synthetic calls may be created"
    pairs = {(c.gene_a, c.gene_b) for c in calls}
    assert ("BCL6", "IGH") not in pairs
    assert ("BCL6", "IGL") not in pairs
    assert ("BCL6", "IGK") not in pairs
    assert not any(c.known_partner for c in calls)


def test_artefact_masked_breakends_are_flagged_and_capped():
    """An unresolvable breakend is reported honestly, never as a resolved fusion."""
    bcl6 = _fc("S1", "3", 187_700_000, "2", 32_916_300, gene_a="BCL6",
               tier="T1", sr=3000, pe=1000)
    bcl6.qc_flags.append("builtin_artefact_locus")
    flag_artefact_masked_breakends([bcl6], cfg=RescueConfig())
    assert "partner_undetermined" in bcl6.qc_flags
    assert bcl6.tier == "T3", "an unresolvable breakend must not hold a T1 tier"
    assert bcl6.gene_b == "", "no partner may be named"


def test_non_artefact_calls_are_untouched():
    """A real measured junction must not be demoted by the artefact handler."""
    real = _fc("S1", "14", 106_500_000, "18", 63_100_000,
               gene_a="IGH", gene_b="BCL2", tier="T1", sr=40, pe=37)
    flag_artefact_masked_breakends([real], cfg=RescueConfig())
    assert real.tier == "T1"
    assert "partner_undetermined" not in real.qc_flags


# ---- the noise filter that starves the artefact channel ----

@pytest.mark.parametrize("seq", [
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTGCAATGAGGTGGGGGGGGGGGGGGGGGG",  # adapter read-through
    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",                       # 2-colour poly-G tail
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",                       # poly-C (reverse strand)
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCACCACCACATTGCACACTCTTTCCCTACAC",  # poly-C + adapter RC
])
def test_noise_clips_are_rejected(seq):
    assert is_noise_clip(seq) is True


@pytest.mark.parametrize("seq", [
    "ACGTACGTTGACCTAGCATCGATCGGATCAAGCTTAACGTACGTACGTTAGCATCAGCATG",  # ordinary genomic
    "TTAGGCATCGATCGAAGCTTGCATCGATCAGCTAGCTAGCATCGATCGATCAGCTAGCTAG",
])
def test_genomic_clips_are_kept(seq):
    assert is_noise_clip(seq) is False


# ---- interchromosomal calls need a read pair across the junction ----

def test_real_translocation_survives_with_pe_and_no_split_reads():
    """The canonical shape of a genuine IG translocation: PE-rich, SR=0.

    9 of 13 confirmed events in the validation cohort look exactly like this —
    the IG switch region is too repetitive to place split reads, but the mates
    straddle the junction and map independently of the junction sequence.
    """
    real = _fc("S1", "14", 106_500_000, "18", 63_100_000,
               gene_a="IGH", gene_b="BCL2", tier="T1", sr=0, pe=21)
    assert demote_interchromosomal_without_discordant([real]) == 0
    assert real.tier == "T1"


def test_split_read_only_interchromosomal_is_demoted():
    """The artefact shape: split reads but no pair across the junction.

    A repeat cross-maps the clipped read while its mate maps normally nearby, so
    no discordant pair forms. This class recurs at identical coordinates across
    unrelated patients and in germline normals.
    """
    art = _fc("S1", "2", 89_790_000, "16", 46_390_722,
              gene_a="IGK", gene_b="", tier="T2", sr=32, pe=0)
    assert demote_interchromosomal_without_discordant([art]) == 1
    assert art.tier == "T3"
    assert "no_discordant_support" in art.qc_flags


def test_canonical_pair_is_not_exempt_from_the_discordant_requirement():
    """A canonical gene pair with no pair across the junction is the fabrication
    pattern itself, so being canonical must not buy an exemption."""
    fake = _fc("S1", "3", 187_700_000, "14", 106_500_000,
               gene_a="BCL6", gene_b="IGH", tier="T1", sr=40, pe=0)
    fake.known_partner = True
    assert demote_interchromosomal_without_discordant([fake]) == 1
    assert fake.tier == "T3"


def test_assembly_supported_calls_are_exempt():
    """An assembled contig across the junction is not a cross-mapped clip.

    Fairness matters here as much as correctness: assembly-based callers report
    PE=0 by construction for their strongest calls (svaba sets dr_evi=0 when
    EVDNC=ASSMB; factera hardcodes discordant_pairs=0). Without this exemption
    the rule would delete a competitor's best evidence and flatter our own tool.
    """
    asm = _fc("S1", "8", 127_740_000, "14", 106_500_000,
              gene_a="MYC", gene_b="IGH", tier="T1", sr=20, pe=0)
    asm.assembly_contigs = 1
    assert demote_interchromosomal_without_discordant([asm]) == 0
    assert asm.tier == "T1"


def test_our_own_scanner_gets_no_assembly_exemption():
    """The scanner always emits assembly_contigs=0, so the rule binds us fully."""
    ours = _fc("S1", "2", 89_790_000, "16", 46_390_722,
               gene_a="IGK", gene_b="", tier="T2", sr=32, pe=0)
    ours.assembly_contigs = 0
    assert demote_interchromosomal_without_discordant([ours]) == 1
    assert ours.tier == "T3"


def test_intrachromosomal_calls_are_untouched():
    """A short-range event sits inside one insert, so no discordant pair is
    expected and its absence is not evidence of an artefact."""
    intra = _fc("S1", "3", 187_700_000, "3", 187_700_400,
                gene_a="BCL6", gene_b="BCL6", tier="T2", sr=12, pe=0)
    assert demote_interchromosomal_without_discordant([intra]) == 0
    assert intra.tier == "T2"


def test_manta_bnd_depth_is_not_an_assembly_signal():
    """BND_DEPTH must not exempt Manta from the discordant-support cull.

    BND_DEPTH is a read-depth annotation Manta puts on EVERY BND record in
    diploidSV (60/60 across the cohort's VCFs; 0/62 in candidateSV). Treating it
    as an assembly signal exempted every Manta breakend from the false-positive
    cull that our own calls all receive — and made the exemption depend on which
    output file the harness parsed rather than on the evidence.
    """
    import tempfile
    from quasarsv.parsers.manta import parse_manta
    vcf = (
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "chr14\t106500000\tA\tN\tN[chr18:63100000[\t99\tPASS\t"
        "SVTYPE=BND;MATEID=B;BND_DEPTH=42\tGT:PR:SR\t0/1:10,4:10,2\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as fh:
        fh.write(vcf)
        p = fh.name
    calls = parse_manta(p, "S")
    assert calls, "record should parse"
    assert calls[0].evidence.assembly_contigs == 0, (
        "BND_DEPTH is depth, not assembly — it must not confer an exemption")


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
