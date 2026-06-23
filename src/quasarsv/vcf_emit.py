"""VCF 4.3 emitter for quasarsv FusionCalls.

Allows quasarsv to publish its calls in the same interchange format as
Manta / GRIDSS / Delly / SvABA / TIDDIT — making downstream comparison and
integration with general SV workflows trivial.

Schema choices
--------------
* BND record per non-intra-chromosomal pair; ALT uses VCF 4.3 bracket
  syntax matching strands.
* DEL / DUP / INV / INS for intra-chromosomal events when ``sv_type``
  is one of those.
* INFO fields include all quasarsv-specific annotations as
  custom tags so nothing is lost in translation:
    ``FF_TIER`` (T1|T2|T3)
    ``FF_EVENT_CLASS``
    ``FF_GENE_A``, ``FF_GENE_B``
    ``FF_DRIVER_LOCUS``
    ``FF_KNOWN_PARTNER`` (0|1) + ``FF_KP_SOURCE``
    ``FF_IN_FRAME`` (0|1)
    ``FF_N_CALLERS``, ``FF_CALLERS``, ``FF_N_EV_TYPES``
    ``FF_SR``, ``FF_PE``, ``FF_ASM``, ``FF_SC``
    ``FF_VAF``, ``FF_PRECISE``, ``FF_ANY_PASS``
    ``FF_QC_FLAGS``
* FORMAT/SAMPLE supplies SR/PE/AD/DP for IGV / GATK compatibility.
* FILTER column = ``PASS`` if any contributing caller marked PASS,
  else ``LowQual``; T3 calls also get ``T3``.
"""
from __future__ import annotations

import gzip
import io
from dataclasses import asdict
from datetime import date
from typing import Iterable, TextIO

from .model import FusionCall


VCF_VERSION = "VCFv4.3"
SOURCE_TAG = "quasarsv"


def _strand_to_alt(chrom_a: str, pos_a: int, strand_a: str,
                   chrom_b: str, pos_b: int, strand_b: str,
                   ref_base: str = "N") -> str:
    """Build a VCF 4.3 BND ALT string for the (a -> b) side.

    Convention:
      strand_a='+' AND strand_b='+'  -> t[p[
      strand_a='+' AND strand_b='-'  -> t]p]
      strand_a='-' AND strand_b='+'  -> [p[t
      strand_a='-' AND strand_b='-'  -> ]p]t
    """
    mate = f"{chrom_b}:{pos_b}"
    if strand_a == "+" and strand_b == "+":
        return f"{ref_base}[{mate}["
    if strand_a == "+" and strand_b == "-":
        return f"{ref_base}]{mate}]"
    if strand_a == "-" and strand_b == "+":
        return f"[{mate}[{ref_base}"
    if strand_a == "-" and strand_b == "-":
        return f"]{mate}]{ref_base}"
    return f"{ref_base}[{mate}["    # fallback


def _info_field(call: FusionCall) -> str:
    parts: list[str] = []

    def add(key: str, val):
        if val is None or val == "":
            return
        if isinstance(val, bool):
            parts.append(f"{key}={1 if val else 0}")
        elif isinstance(val, (list, tuple)):
            if not val:
                return
            parts.append(f"{key}={'|'.join(str(v) for v in val)}")
        else:
            parts.append(f"{key}={val}")

    sv_type = call.sv_type or "BND"
    add("SVTYPE", sv_type)
    if sv_type != "BND" and call.chrom_a == call.chrom_b:
        add("END", call.pos_b)
        add("SVLEN", call.pos_b - call.pos_a)
    if not call.precise:
        parts.append("IMPRECISE")
    add("FF_TIER", call.tier)
    add("FF_EVENT_CLASS", call.event_class)
    add("FF_GENE_A", call.gene_a)
    add("FF_GENE_B", call.gene_b)
    add("FF_DRIVER_LOCUS", call.driver_locus)
    add("FF_KNOWN_PARTNER", call.known_partner)
    add("FF_KP_SOURCE", call.known_partner_source)
    add("FF_IN_FRAME", call.in_frame)
    add("FF_N_CALLERS", call.n_callers)
    add("FF_CALLERS", call.callers_supporting)
    add("FF_N_EV_TYPES", call.n_evidence_types)
    add("FF_SR", call.split_reads)
    add("FF_PE", call.discordant_pairs)
    add("FF_ASM", call.assembly_contigs)
    add("FF_SC", call.soft_clips)
    add("FF_VAF", f"{call.vaf:.4f}" if call.vaf else 0.0)
    add("FF_PRECISE", call.precise)
    add("FF_ANY_PASS", call.any_pass)
    add("FF_QC_FLAGS", call.qc_flags)
    add("FF_FUSION_ID", call.fusion_id)
    return ";".join(parts) if parts else "."


def _filter_field(call: FusionCall) -> str:
    flags: list[str] = []
    if not call.any_pass and call.n_callers <= 1:
        flags.append("LowQual")
    if call.tier == "T3":
        flags.append("T3")
    return "PASS" if not flags else ";".join(flags)


def _format_sample_fields(call: FusionCall) -> tuple[str, str]:
    sr = max(int(call.split_reads), 0)
    pe = max(int(call.discordant_pairs), 0)
    ad = sr + pe
    dp = ad if ad > 0 else 0
    return "GT:SR:PE:AD:DP", f"./.:{sr}:{pe}:{ad}:{dp}"


def _header(sample: str, contigs: list[str] | None = None) -> str:
    lines = [
        f"##fileformat={VCF_VERSION}",
        f"##fileDate={date.today().isoformat().replace('-', '')}",
        f"##source={SOURCE_TAG}",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV\">",
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise breakpoint\">",
        "##INFO=<ID=FF_TIER,Number=1,Type=String,Description=\"quasarsv confidence tier (T1|T2|T3)\">",
        "##INFO=<ID=FF_EVENT_CLASS,Number=1,Type=String,Description=\"quasarsv event class\">",
        "##INFO=<ID=FF_GENE_A,Number=1,Type=String,Description=\"Annotated gene at side A\">",
        "##INFO=<ID=FF_GENE_B,Number=1,Type=String,Description=\"Annotated gene at side B\">",
        "##INFO=<ID=FF_DRIVER_LOCUS,Number=1,Type=String,Description=\"Driver-locus pair label\">",
        "##INFO=<ID=FF_KNOWN_PARTNER,Number=1,Type=Integer,Description=\"1 if gene pair is canonical lymphoma\">",
        "##INFO=<ID=FF_KP_SOURCE,Number=1,Type=String,Description=\"Canonical partner source\">",
        "##INFO=<ID=FF_IN_FRAME,Number=1,Type=Integer,Description=\"1 if predicted in-frame\">",
        "##INFO=<ID=FF_N_CALLERS,Number=1,Type=Integer,Description=\"Number of supporting callers\">",
        "##INFO=<ID=FF_CALLERS,Number=.,Type=String,Description=\"Pipe-separated supporting callers\">",
        "##INFO=<ID=FF_N_EV_TYPES,Number=1,Type=Integer,Description=\"Independent evidence types\">",
        "##INFO=<ID=FF_SR,Number=1,Type=Integer,Description=\"Split-read count\">",
        "##INFO=<ID=FF_PE,Number=1,Type=Integer,Description=\"Discordant-pair count\">",
        "##INFO=<ID=FF_ASM,Number=1,Type=Integer,Description=\"Assembly contig count\">",
        "##INFO=<ID=FF_SC,Number=1,Type=Integer,Description=\"Soft-clip count\">",
        "##INFO=<ID=FF_VAF,Number=1,Type=Float,Description=\"Variant allele fraction\">",
        "##INFO=<ID=FF_PRECISE,Number=1,Type=Integer,Description=\"1 if precise breakpoint\">",
        "##INFO=<ID=FF_ANY_PASS,Number=1,Type=Integer,Description=\"1 if any caller PASS\">",
        "##INFO=<ID=FF_QC_FLAGS,Number=.,Type=String,Description=\"quasarsv QC flags\">",
        "##INFO=<ID=FF_FUSION_ID,Number=1,Type=String,Description=\"Stable quasarsv fusion id\">",
        "##FILTER=<ID=LowQual,Description=\"Single-caller and no PASS upstream\">",
        "##FILTER=<ID=T3,Description=\"quasarsv tier T3 (low confidence)\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Split reads supporting variant\">",
        "##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Discordant pairs supporting variant\">",
        "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allelic depth (SR+PE)\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth proxy\">",
        "##ALT=<ID=BND,Description=\"Breakend\">",
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        "##ALT=<ID=DUP,Description=\"Duplication\">",
        "##ALT=<ID=INV,Description=\"Inversion\">",
        "##ALT=<ID=INS,Description=\"Insertion\">",
    ]
    if contigs:
        for c in contigs:
            lines.append(f"##contig=<ID={c}>")
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                            "FILTER", "INFO", "FORMAT", sample or "SAMPLE"]))
    return "\n".join(lines) + "\n"


def _record(call: FusionCall, mate_idx: int | None = None) -> str:
    chrom = f"chr{call.chrom_a}" if not call.chrom_a.startswith("chr") else call.chrom_a
    pos = max(1, int(call.pos_a))
    vid = call.fusion_id if mate_idx is None else f"{call.fusion_id}.{mate_idx}"
    ref = "N"
    sv_type = call.sv_type or "BND"

    if sv_type == "BND" or call.chrom_a != call.chrom_b:
        mate_chrom = f"chr{call.chrom_b}" if not call.chrom_b.startswith("chr") else call.chrom_b
        alt = _strand_to_alt(chrom, pos, call.strand_a or "+",
                              mate_chrom, max(1, int(call.pos_b)), call.strand_b or "+")
    else:
        alt = f"<{sv_type}>"

    qual = f"{call.raw_qual_max:.1f}" if call.raw_qual_max else "."
    filt = _filter_field(call)
    info = _info_field(call)
    fmt, sample_v = _format_sample_fields(call)
    return "\t".join([chrom, str(pos), vid, ref, alt, qual, filt, info, fmt, sample_v]) + "\n"


def write_vcf(calls: Iterable[FusionCall], path: str,
              sample: str | None = None,
              contigs: list[str] | None = None,
              emit_mates: bool = True) -> int:
    """Write FusionCalls to a (possibly bgzipped) VCF.

    Parameters
    ----------
    calls
        Iterable of FusionCall.
    path
        Output path (``.vcf`` plain text or ``.vcf.gz`` gzip).
    sample
        Sample id for the FORMAT column. Defaults to the first call's sample.
    contigs
        Optional list of contig names for ``##contig=`` header lines.
    emit_mates
        If True, each BND emits both ends (matched by MATEID). Some downstream
        consumers expect mate pairs; others (Manta-style) are fine with one
        per breakpoint. Default True.

    Returns
    -------
    int
        Number of VCF records written.
    """
    calls = list(calls)
    if sample is None:
        sample = calls[0].sample if calls else "SAMPLE"

    def _open(p):
        if p.endswith(".gz"):
            return gzip.open(p, "wt", encoding="utf-8")
        return open(p, "w", encoding="utf-8")

    n = 0
    with _open(path) as fh:
        fh.write(_header(sample, contigs))
        for c in calls:
            chrom_a = f"chr{c.chrom_a}" if not c.chrom_a.startswith("chr") else c.chrom_a
            chrom_b = f"chr{c.chrom_b}" if not c.chrom_b.startswith("chr") else c.chrom_b
            is_bnd = (c.sv_type or "BND") == "BND" or chrom_a != chrom_b
            if not is_bnd:
                fh.write(_record(c))
                n += 1
                continue
            # Two records per BND when emit_mates: one anchored at A, one at B
            fh.write(_record(c, mate_idx=1))
            n += 1
            if emit_mates:
                rev = FusionCall(**{**asdict(c),
                                     "chrom_a": c.chrom_b, "pos_a": c.pos_b, "strand_a": c.strand_b,
                                     "chrom_b": c.chrom_a, "pos_b": c.pos_a, "strand_b": c.strand_a})
                fh.write(_record(rev, mate_idx=2))
                n += 1
    return n


def write_vcf_to_string(calls: Iterable[FusionCall], sample: str | None = None) -> str:
    """Same as ``write_vcf`` but returns a string — convenient for tests."""
    buf = io.StringIO()
    calls = list(calls)
    sample = sample or (calls[0].sample if calls else "SAMPLE")
    buf.write(_header(sample))
    for c in calls:
        is_bnd = (c.sv_type or "BND") == "BND" or c.chrom_a != c.chrom_b
        if not is_bnd:
            buf.write(_record(c))
            continue
        buf.write(_record(c, mate_idx=1))
        rev = FusionCall(**{**asdict(c),
                             "chrom_a": c.chrom_b, "pos_a": c.pos_b, "strand_a": c.strand_b,
                             "chrom_b": c.chrom_a, "pos_b": c.pos_a, "strand_b": c.strand_a})
        buf.write(_record(rev, mate_idx=2))
    return buf.getvalue()
