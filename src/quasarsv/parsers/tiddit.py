"""TIDDIT VCF parser.

TIDDIT 3.x emits INV, DEL, DUP, INS, TDUP, BND with INFO fields:
  SVTYPE, END (or for BND parsed from ALT), CIPOS, CIEND, REGIONA, REGIONB,
  LFA / LTE / LFB / LTB (per-side discordant / split-read counts).

Quality is in QUAL (per-record). FILTER is PASS / various tiddit-specific tags.
"""
from __future__ import annotations

from ..model import BreakpointCall, Evidence
from .base import iter_records, parse_info, parse_bnd_alt, normalise_order


def _int(v, default=0):
    if v is None or v == "" or v == ".":
        return default
    try:
        return int(v)
    except ValueError:
        return default


def parse_tiddit(path: str, sample: str) -> list[BreakpointCall]:
    seen_mates: set[str] = set()
    out: list[BreakpointCall] = []
    for f in iter_records(path):
        chrom, pos, vid, ref, alt, qual, filt, info_str = f[:8]
        info = parse_info(info_str)
        sv_type = info.get("SVTYPE", "")
        if not sv_type:
            continue
        precise = "IMPRECISE" not in info
        pass_filter = filt in ("PASS", ".")
        try:
            q = float(qual) if qual not in ("", ".") else 0.0
        except ValueError:
            q = 0.0

        # TIDDIT discordant + split counts (LFA = forward-strand discordants on A side,
        # LTE = trans-evidence, etc.). Sum read-pair-style fields conservatively.
        disc = (_int(info.get("LFA")) + _int(info.get("LFB"))
                + _int(info.get("LTA")) + _int(info.get("LTB"))
                + _int(info.get("DV")))
        split = _int(info.get("SR")) + _int(info.get("RV"))

        if sv_type == "BND":
            mate_id = info.get("MATEID", "")
            if vid in seen_mates or (mate_id and mate_id in seen_mates):
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
            end = _int(info.get("END"), int(pos))
            sa = "+"
            sb = "-" if sv_type in ("DEL", "TDUP", "DUP") else ("+" if sv_type == "INV" else "+")
            chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
                chrom, int(pos), sa, chrom, end, sb
            )

        ev = Evidence(
            caller="tiddit",
            split_reads=split,
            discordant_pairs=disc,
            assembly_contigs=0,
            mapq=0,
            vaf=0.0,
            filter_pass=pass_filter,
            precise=precise,
            raw_qual=q,
        )
        out.append(BreakpointCall(
            sample=sample, chrom_a=chrom_a, pos_a=pos_a, strand_a=strand_a,
            chrom_b=chrom_b, pos_b=pos_b, strand_b=strand_b,
            sv_type=sv_type, evidence=ev, record_id=vid,
        ))
    return out
