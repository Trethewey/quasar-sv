"""Delly VCF parser.

Delly emits per-event records (DEL/DUP/INV/INS/BND with CT connection types).
Evidence in INFO: PE (discordant pairs), SR (split reads if assembled);
in FORMAT per-sample: DR/DV (ref/alt discordant), RR/RV (ref/alt split-read).
"""
from __future__ import annotations

from ..model import BreakpointCall, Evidence
from .base import iter_records, parse_info, parse_format_sample, parse_bnd_alt, normalise_order


# Delly CT (connection type) -> strands
CT_STRANDS = {
    "3to5": ("+", "-"),
    "5to3": ("-", "+"),
    "3to3": ("+", "+"),
    "5to5": ("-", "-"),
}


def parse_delly(path: str, sample: str) -> list[BreakpointCall]:
    out: list[BreakpointCall] = []
    for f in iter_records(path):
        chrom, pos, vid, ref, alt, qual, filt, info_str = f[:8]
        info = parse_info(info_str)
        sv_type = info.get("SVTYPE", "")
        if not sv_type:
            continue
        precise = "IMPRECISE" not in info
        pass_filter = (filt == "PASS")
        try:
            q = float(qual) if qual not in ("", ".") else 0.0
        except ValueError:
            q = 0.0

        # Default strands by CT
        ct = info.get("CT", "3to5")
        sa, sb = CT_STRANDS.get(ct, ("+", "-"))

        if sv_type == "BND":
            parsed = parse_bnd_alt(alt)
            if parsed is None:
                continue
            mchrom, mpos, sa, sb = parsed
            chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
                chrom, int(pos), sa, mchrom, mpos, sb
            )
        else:
            end = int(info.get("END", pos))
            chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
                chrom, int(pos), sa, chrom, end, sb
            )

        # Evidence
        pe_info = int(info.get("PE") or 0)
        sr_info = int(info.get("SR") or 0)
        fmt_sample = {}
        if len(f) >= 10:
            fmt_sample = parse_format_sample(f[8], f[9])
        dv = int(fmt_sample.get("DV") or 0)
        rv = int(fmt_sample.get("RV") or 0)
        dr = int(fmt_sample.get("DR") or 0)
        rr = int(fmt_sample.get("RR") or 0)
        depth = dv + dr + rv + rr
        vaf = ((dv + rv) / depth) if depth > 0 else 0.0
        mq = int(info.get("MAPQ") or 0)
        ci = _parse_ci(info.get("CIPOS"))
        ci_end = _parse_ci(info.get("CIEND"))

        ev = Evidence(
            caller="delly",
            split_reads=max(sr_info, rv),
            discordant_pairs=max(pe_info, dv),
            assembly_contigs=int(info.get("CONSBP") is not None or info.get("CONSENSUS") is not None),
            mapq=mq,
            vaf=vaf,
            filter_pass=pass_filter,
            precise=precise,
            raw_qual=q,
        )
        out.append(BreakpointCall(
            sample=sample, chrom_a=chrom_a, pos_a=pos_a, strand_a=strand_a,
            chrom_b=chrom_b, pos_b=pos_b, strand_b=strand_b,
            sv_type=sv_type, evidence=ev,
            cipos_a=ci, cipos_b=ci_end, record_id=vid,
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
