"""Tests for the VCF emitter."""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from quasarsv.model import FusionCall
from quasarsv.vcf_emit import write_vcf, write_vcf_to_string


def _fc(sample="S1", chrom_a="14", pos_a=106_000_000, strand_a="+",
        chrom_b="18", pos_b=63_200_000, strand_b="-",
        sv_type="BND", gene_a="IGH", gene_b="BCL2",
        tier="T1", known_partner=True, kp_src="canonical:Follicular",
        split_reads=12, discordant_pairs=24, qc_flags=None):
    return FusionCall(
        sample=sample,
        fusion_id=f"{sample}__{gene_a}_{gene_b}",
        chrom_a=chrom_a, pos_a=pos_a, strand_a=strand_a,
        chrom_b=chrom_b, pos_b=pos_b, strand_b=strand_b,
        sv_type=sv_type,
        gene_a=gene_a, gene_b=gene_b,
        tier=tier, known_partner=known_partner, known_partner_source=kp_src,
        split_reads=split_reads, discordant_pairs=discordant_pairs,
        callers_supporting=["manta", "quasar"], n_callers=2,
        n_evidence_types=2, precise=True, any_pass=True,
        qc_flags=qc_flags or [],
    )


def test_vcf_emits_header_and_record():
    text = write_vcf_to_string([_fc()])
    assert text.startswith("##fileformat=VCFv4.3")
    assert "##INFO=<ID=FF_TIER" in text
    assert "##FORMAT=<ID=SR" in text
    body = [l for l in text.splitlines() if not l.startswith("#")]
    assert any("FF_TIER=T1" in line for line in body)
    assert any("FF_KNOWN_PARTNER=1" in line for line in body)


def test_vcf_bnd_alt_syntax():
    text = write_vcf_to_string([
        _fc(strand_a="+", strand_b="+", chrom_b="18", pos_b=63_200_000),
    ])
    # strand_a='+', strand_b='+' -> N[chr18:63200000[
    assert "N[chr18:63200000[" in text


def test_vcf_bnd_emits_mate_pair_by_default():
    text = write_vcf_to_string([_fc()])
    body = [l for l in text.splitlines() if l and not l.startswith("#")]
    # Two records expected (mate + reverse)
    assert len(body) == 2


def test_vcf_filter_field_uses_passes_and_tier():
    pass_call = _fc(tier="T1")
    t3_call = _fc(tier="T3", known_partner=False, kp_src="", split_reads=2, discordant_pairs=0)
    t3_call.any_pass = False
    t3_call.n_callers = 1
    txt_pass = write_vcf_to_string([pass_call])
    txt_t3 = write_vcf_to_string([t3_call])
    body_pass = [l for l in txt_pass.splitlines() if l and not l.startswith("#")][0].split("\t")
    body_t3 = [l for l in txt_t3.splitlines() if l and not l.startswith("#")][0].split("\t")
    assert body_pass[6] == "PASS"     # FILTER column
    assert "T3" in body_t3[6]


def test_vcf_round_trips_through_file():
    call = _fc()
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.vcf"
        n = write_vcf([call], str(out))
        assert n == 2     # mate pair
        text = out.read_text(encoding="utf-8")
        assert "##fileformat=VCFv4.3" in text
        assert "FF_GENE_A=IGH" in text
        assert "FF_GENE_B=BCL2" in text


def test_vcf_gz_extension_is_gzipped():
    call = _fc()
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.vcf.gz"
        write_vcf([call], str(out))
        import gzip
        with gzip.open(out, "rt") as fh:
            text = fh.read()
        assert "##fileformat=VCFv4.3" in text


def test_vcf_intra_chr_emits_sv_type_alt():
    """Intra-chromosomal DEL gets <DEL> ALT, not BND."""
    call = _fc(chrom_a="3", chrom_b="3", pos_a=187_500_000, pos_b=187_700_000,
               sv_type="DEL")
    text = write_vcf_to_string([call])
    body = [l for l in text.splitlines() if l and not l.startswith("#")]
    assert len(body) == 1
    fields = body[0].split("\t")
    assert fields[4] == "<DEL>"


def test_vcf_qc_flags_preserved():
    call = _fc(qc_flags=["inferred_via_artefact_rescue", "ig_partner_ambiguous"])
    text = write_vcf_to_string([call])
    assert "FF_QC_FLAGS=inferred_via_artefact_rescue|ig_partner_ambiguous" in text


def _body(text):
    return [l.split("\t") for l in text.splitlines() if l and not l.startswith("#")]


def test_vcf_bnd_mate_reciprocal_orientation():
    # VCF 4.3 §5.4: t[p[  <->  ]p]t
    text = write_vcf_to_string([_fc(strand_a="+", strand_b="+",
                                    chrom_a="14", pos_a=106_000_000,
                                    chrom_b="18", pos_b=63_200_000)])
    alt = {r[0]: r[4] for r in _body(text)}
    assert alt["chr14"] == "N[chr18:63200000["
    assert alt["chr18"] == "]chr14:106000000]N"


def test_vcf_all_four_strand_combos_round_trip():
    from quasarsv.parsers.base import parse_bnd_alt
    for sa, sb in [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")]:
        text = write_vcf_to_string([_fc(strand_a=sa, strand_b=sb,
                                        chrom_a="14", pos_a=106_000_000,
                                        chrom_b="18", pos_b=63_200_000)])
        by = {r[0]: r for r in _body(text)}
        prim = parse_bnd_alt(by["chr14"][4])
        assert prim[:2] == ("chr18", 63_200_000) and (prim[2], prim[3]) == (sa, sb)
        mate = parse_bnd_alt(by["chr18"][4])
        assert mate[:2] == ("chr14", 106_000_000)   # mate points back at the primary


def test_vcf_mateid_cross_references():
    recs = _body(write_vcf_to_string([_fc()]))
    ids = {r[2] for r in recs}
    assert len(ids) == 2
    for r in recs:
        info = dict(kv.split("=", 1) for kv in r[7].split(";") if "=" in kv)
        assert info["MATEID"] in ids and info["MATEID"] != r[2]


def test_vcf_records_coordinate_sorted():
    c1 = _fc(chrom_a="18", pos_a=63_000_000, chrom_b="2", pos_b=30_000_000)
    c2 = _fc(chrom_a="2", pos_a=10_000_000, chrom_b="8", pos_b=127_000_000)
    recs = _body(write_vcf_to_string([c1, c2]))

    def ck(c):
        s = c[3:]
        return int(s) if s.isdigit() else 99
    keys = [(ck(r[0]), int(r[1])) for r in recs]
    assert keys == sorted(keys)


def test_vcf_info_values_percent_encoded():
    text = write_vcf_to_string([_fc(kp_src="canonical:Burkitt; DH-DLBCL")])
    info = _body(text)[0][7]
    assert " " not in info                       # no raw whitespace in INFO
    assert "Burkitt%3B%20DH-DLBCL" in info        # ';' and ' ' encoded


def test_vcf_svlen_negative_for_del():
    call = _fc(chrom_a="3", chrom_b="3", pos_a=187_500_000, pos_b=187_700_000, sv_type="DEL")
    info = _body(write_vcf_to_string([call]))[0][7]
    assert "SVLEN=-200000" in info


def test_vcf_has_contig_headers():
    text = write_vcf_to_string([_fc()])
    assert "##contig=<ID=chr14>" in text
    assert "##contig=<ID=chr18>" in text
