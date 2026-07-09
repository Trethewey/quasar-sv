"""VCF 4.3 emitter for quasarsv FusionCalls.

Allows quasarsv to publish its calls in the same interchange format as
Manta / GRIDSS / Delly / SvABA / TIDDIT — making downstream comparison and
integration with general SV workflows trivial.

Spec compliance
---------------
* One BND record per breakend; a translocation emits two mated records whose
  ALTs form a valid VCF 4.3 reciprocal pair, linked by ``MATEID``.
* Records are written in coordinate order (contig, then POS) so the output is
  ``tabix``/``bcftools`` indexable; ``.vcf.gz`` output is **bgzip**-compressed
  when pysam is available (falls back to gzip otherwise).
* INFO string values are percent-encoded (``;`` ``=`` ``,`` whitespace ``%``),
  so free-text like a partner-source disease label cannot corrupt the column.
* ``SVLEN`` is negative for DEL, per the VCF 4.3 sign convention.

Schema choices
--------------
* BND record per non-intra-chromosomal pair; ALT uses VCF 4.3 bracket syntax.
* DEL / DUP / INV / INS for intra-chromosomal events when ``sv_type`` is one.
* INFO carries all quasarsv-specific annotation as ``FF_*`` custom tags.
* FORMAT/SAMPLE supplies SR/PE/AD/DP for IGV / GATK compatibility.
* FILTER = ``PASS`` if any contributing caller passed, else ``LowQual``;
  T3 calls also get ``T3``.
"""
from __future__ import annotations

import gzip
import io
from dataclasses import asdict
from datetime import date, timezone, datetime
from typing import Iterable

from .model import FusionCall


VCF_VERSION = "VCFv4.3"
SOURCE_TAG = "quasarsv"

# VCF 4.3 reciprocal-breakend rule, expressed in this module's (strand_a,
# strand_b) encoding (see _strand_to_alt). The mate of a breakend uses these
# strands, points back at the primary's coordinates, and carries the primary's
# reference base. Verified against the VCF 4.3 §5.4 canonical example:
#   t[p[ <-> ]p]t ,  ]p]t <-> t[p[ ,  t]p] <-> t]p] ,  [p[t <-> [p[t
_MATE_STRANDS = {
    ("+", "+"): ("-", "-"),
    ("-", "-"): ("+", "+"),
    ("+", "-"): ("+", "-"),
    ("-", "+"): ("-", "+"),
}


def _pct(value) -> str:
    """Percent-encode INFO-reserved characters in a string value (VCF 4.3)."""
    s = str(value)
    for ch, rep in (("%", "%25"), (";", "%3B"), ("=", "%3D"), (",", "%2C"),
                    (" ", "%20"), ("\t", "%09"), ("\n", "%0A"), ("\r", "%0D")):
        s = s.replace(ch, rep)
    return s


def _disp_chrom(chrom: str) -> str:
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def _contig_sort_key(chrom: str) -> tuple[int, int, str]:
    """Order contigs 1..22, X, Y, then everything else alphabetically."""
    s = chrom[3:] if chrom.startswith("chr") else chrom
    try:
        return (0, int(s), "")
    except ValueError:
        special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        if s in special:
            return (0, special[s], "")
        return (1, 0, s)


def _strand_to_alt(mate_chrom: str, mate_pos: int,
                   strand_a: str, strand_b: str, ref_base: str = "N") -> str:
    """Build a VCF 4.3 BND ALT for a breakend joining to ``mate_chrom:mate_pos``.

    Convention (inverse of parsers.base.parse_bnd_alt):
      strand_a='+' AND strand_b='+'  -> t[p[
      strand_a='+' AND strand_b='-'  -> t]p]
      strand_a='-' AND strand_b='+'  -> [p[t
      strand_a='-' AND strand_b='-'  -> ]p]t
    """
    mate = f"{mate_chrom}:{mate_pos}"
    if strand_a == "+" and strand_b == "+":
        return f"{ref_base}[{mate}["
    if strand_a == "+" and strand_b == "-":
        return f"{ref_base}]{mate}]"
    if strand_a == "-" and strand_b == "+":
        return f"[{mate}[{ref_base}"
    if strand_a == "-" and strand_b == "-":
        return f"]{mate}]{ref_base}"
    return f"{ref_base}[{mate}["    # fallback


def _info_field(call: FusionCall, *, swap_genes: bool = False,
                mate_id: str | None = None) -> str:
    """Build the INFO column. ``swap_genes`` puts the local breakend's gene in
    FF_GENE_A (used for the mate record, whose local side is B)."""
    parts: list[str] = []

    def add(key: str, val, *, encode: bool = False):
        if val is None or val == "":
            return
        if isinstance(val, bool):
            parts.append(f"{key}={1 if val else 0}")
        elif isinstance(val, (list, tuple)):
            if not val:
                return
            parts.append(f"{key}={'|'.join(_pct(v) for v in val)}")
        else:
            parts.append(f"{key}={_pct(val) if encode else val}")

    gene_a, gene_b = (call.gene_b, call.gene_a) if swap_genes else (call.gene_a, call.gene_b)

    sv_type = call.sv_type or "BND"
    add("SVTYPE", sv_type)
    if mate_id:
        add("MATEID", mate_id)
    if sv_type != "BND" and call.chrom_a == call.chrom_b:
        add("END", call.pos_b)
        span = abs(int(call.pos_b) - int(call.pos_a))
        add("SVLEN", -span if sv_type == "DEL" else span)
    if not call.precise:
        parts.append("IMPRECISE")
    add("FF_TIER", call.tier)
    add("FF_EVENT_CLASS", call.event_class, encode=True)
    add("FF_GENE_A", gene_a, encode=True)
    add("FF_GENE_B", gene_b, encode=True)
    add("FF_DRIVER_LOCUS", call.driver_locus, encode=True)
    add("FF_KNOWN_PARTNER", call.known_partner)
    add("FF_KP_SOURCE", call.known_partner_source, encode=True)
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
    add("FF_FUSION_ID", call.fusion_id, encode=True)
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
    # AD is reserved (Number=R): ref support unknown ('.'), alt support = ad.
    return "GT:SR:PE:AD:DP", f"./.:{sr}:{pe}:.,{ad}:{dp}"


def _header(sample: str, contigs: list[str] | None = None) -> str:
    lines = [
        f"##fileformat={VCF_VERSION}",
        f"##fileDate={date.today().isoformat().replace('-', '')}",
        f"##source={SOURCE_TAG}",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">",
        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of mate breakend\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT (negative for DEL)\">",
        "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise breakpoint\">",
        "##INFO=<ID=FF_TIER,Number=1,Type=String,Description=\"quasarsv confidence tier (T1|T2|T3)\">",
        "##INFO=<ID=FF_EVENT_CLASS,Number=1,Type=String,Description=\"quasarsv event class\">",
        "##INFO=<ID=FF_GENE_A,Number=1,Type=String,Description=\"Annotated gene at the local breakend\">",
        "##INFO=<ID=FF_GENE_B,Number=1,Type=String,Description=\"Annotated gene at the mate breakend\">",
        "##INFO=<ID=FF_DRIVER_LOCUS,Number=1,Type=String,Description=\"Driver-locus pair label\">",
        "##INFO=<ID=FF_KNOWN_PARTNER,Number=1,Type=Integer,Description=\"1 if gene pair is canonical lymphoma\">",
        "##INFO=<ID=FF_KP_SOURCE,Number=1,Type=String,Description=\"Canonical partner source (percent-encoded)\">",
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
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for ref and alt (ref unknown)\">",
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


def _records_for_call(call: FusionCall, emit_mates: bool) -> list[tuple]:
    """Return one or two record tuples (sort-ready) for a call.

    Tuple = (contig_key, pos, chrom_display, pos, id, ref, alt, qual, filter,
    info, format, sample_value).
    """
    ref = "N"
    qual = f"{call.raw_qual_max:.1f}" if call.raw_qual_max else "."
    filt = _filter_field(call)
    fmt, sample_v = _format_sample_fields(call)
    sv_type = call.sv_type or "BND"

    ca, cb = _disp_chrom(call.chrom_a), _disp_chrom(call.chrom_b)
    pa, pb = max(1, int(call.pos_a)), max(1, int(call.pos_b))
    sa, sb = call.strand_a or "+", call.strand_b or "+"

    is_bnd = sv_type == "BND" or call.chrom_a != call.chrom_b

    def _tuple(chrom, pos, vid, alt, info):
        return (_contig_sort_key(chrom), pos, vid,
                (chrom, str(pos), vid, ref, alt, qual, filt, info, fmt, sample_v))

    if not is_bnd:
        info = _info_field(call, mate_id=None)
        return [_tuple(ca, pa, call.fusion_id, f"<{sv_type}>", info)]

    # Breakend: primary at A pointing to B.
    if not emit_mates:
        alt = _strand_to_alt(cb, pb, sa, sb, ref)
        info = _info_field(call, mate_id=None)
        return [_tuple(ca, pa, call.fusion_id, alt, info)]

    id1, id2 = f"{call.fusion_id}.1", f"{call.fusion_id}.2"
    alt1 = _strand_to_alt(cb, pb, sa, sb, ref)
    info1 = _info_field(call, swap_genes=False, mate_id=id2)
    # Mate at B pointing to A, with the reciprocal bracket orientation.
    msa, msb = _MATE_STRANDS.get((sa, sb), (sb, sa))
    alt2 = _strand_to_alt(ca, pa, msa, msb, ref)
    info2 = _info_field(call, swap_genes=True, mate_id=id1)
    return [_tuple(ca, pa, id1, alt1, info1),
            _tuple(cb, pb, id2, alt2, info2)]


def _build_sorted_lines(calls: list[FusionCall], emit_mates: bool) -> tuple[list[str], list[str]]:
    """Return (contig_ids, record_lines) with records coordinate-sorted."""
    recs: list[tuple] = []
    for c in calls:
        recs.extend(_records_for_call(c, emit_mates))
    recs.sort(key=lambda r: (r[0], r[1], r[2]))
    contigs = sorted({r[3][0] for r in recs}, key=_contig_sort_key)
    lines = ["\t".join(r[3]) for r in recs]
    return contigs, lines


def write_vcf(calls: Iterable[FusionCall], path: str,
              sample: str | None = None,
              contigs: list[str] | None = None,
              emit_mates: bool = True) -> int:
    """Write FusionCalls to a coordinate-sorted, indexable VCF.

    ``.vcf.gz`` output is bgzip-compressed when pysam is available (so
    ``tabix -p vcf`` / ``bcftools index`` work); otherwise it falls back to
    gzip. Returns the number of VCF records written.
    """
    calls = list(calls)
    if sample is None:
        sample = calls[0].sample if calls else "SAMPLE"

    derived_contigs, lines = _build_sorted_lines(calls, emit_mates)
    text = _header(sample, contigs or derived_contigs or None) + \
        ("".join(l + "\n" for l in lines))

    if path.endswith(".gz"):
        data = text.encode("utf-8")
        try:
            import pysam  # bgzip: gzip-compatible AND tabix-indexable
            with pysam.BGZFile(path, "wb") as fh:
                fh.write(data)
        except Exception:
            with gzip.open(path, "wb") as fh:
                fh.write(data)
    else:
        with open(path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
    return len(lines)


def write_vcf_to_string(calls: Iterable[FusionCall], sample: str | None = None,
                        emit_mates: bool = True) -> str:
    """Same as ``write_vcf`` but returns the (coordinate-sorted) VCF as a string."""
    calls = list(calls)
    sample = sample or (calls[0].sample if calls else "SAMPLE")
    derived_contigs, lines = _build_sorted_lines(calls, emit_mates)
    buf = io.StringIO()
    buf.write(_header(sample, derived_contigs or None))
    for l in lines:
        buf.write(l + "\n")
    return buf.getvalue()
