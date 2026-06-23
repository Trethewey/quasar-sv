"""Manta VCF parser.

Manta emits BND mate pairs (one record per breakpoint end) + DEL/DUP/INV/INS
single-record events. FORMAT carries PR (discordant pairs) and SR (split reads)
as comma-separated ref,alt counts.
"""
from __future__ import annotations

from ..model import BreakpointCall, Evidence
from .base import iter_records, parse_info, parse_format_sample, parse_bnd_alt, normalise_order


def _alt_counts(field: str | None) -> int:
    if not field or field == ".":
        return 0
    parts = field.split(",")
    return int(parts[-1]) if parts[-1].isdigit() else 0


def parse_manta(path: str, sample: str) -> list[BreakpointCall]:
    seen_mates: set[str] = set()
    out: list[BreakpointCall] = []

    for f in iter_records(path):
        chrom, pos, vid, ref, alt, qual, filt, info_str = f[:8]
        info = parse_info(info_str)
        sv_type = info.get("SVTYPE", "")
        if not sv_type:
            continue   # not an SV record (e.g. MuTect2 SNV input)
        precise = "IMPRECISE" not in info
        pass_filter = filt in ("PASS", ".")

        # FORMAT / sample (Manta tumour-only has one sample column)
        fmt_sample = {}
        if len(f) >= 10:
            fmt_sample = parse_format_sample(f[8], f[9])

        split_reads = _alt_counts(fmt_sample.get("SR"))
        disc_pairs = _alt_counts(fmt_sample.get("PR"))

        try:
            q = float(qual) if qual not in ("", ".") else 0.0
        except ValueError:
            q = 0.0

        if sv_type == "BND":
            mate_id = info.get("MATEID", "")
            if vid in seen_mates or mate_id in seen_mates:
                continue
            seen_mates.add(vid)
            parsed = parse_bnd_alt(alt)
            if parsed is None:
                continue
            mchrom, mpos, sa, sb = parsed
            chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
                chrom, int(pos), sa, mchrom, mpos, sb
            )
        else:
            end = int(info.get("END", pos))
            sa = "+"
            sb = "-" if sv_type in ("DEL", "DUP") else ("+" if sv_type == "INV" else "+")
            chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
                chrom, int(pos), sa, chrom, end, sb
            )

        ci_a = _parse_ci(info.get("CIPOS"))
        ci_b = _parse_ci(info.get("CIEND") or info.get("CIPOS"))

        ev = Evidence(
            caller="manta",
            split_reads=split_reads,
            discordant_pairs=disc_pairs,
            assembly_contigs=1 if "CONTIG" in info or "BND_DEPTH" in info else 0,
            mapq=int(info.get("MAPQ", 0) or 0),
            vaf=_safe_float(info.get("BND_DEPTH"), 0.0),
            filter_pass=pass_filter,
            precise=precise,
            raw_qual=q,
        )
        out.append(BreakpointCall(
            sample=sample, chrom_a=chrom_a, pos_a=pos_a, strand_a=strand_a,
            chrom_b=chrom_b, pos_b=pos_b, strand_b=strand_b,
            sv_type=sv_type, evidence=ev,
            cipos_a=ci_a, cipos_b=ci_b, record_id=vid,
        ))
    return out


def _parse_ci(s: str | None) -> tuple[int, int]:
    if not s:
        return (0, 0)
    parts = s.split(",")
    if len(parts) != 2:
        return (0, 0)
    try:
        return (int(parts[0]), int(parts[1]))
    except ValueError:
        return (0, 0)


def _safe_float(s: str | None, default: float) -> float:
    if s is None or s == "":
        return default
    try:
        return float(s)
    except ValueError:
        return default
