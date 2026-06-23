"""SvABA VCF parser.

SvABA outputs are paired BND records (MATEID). INFO carries EVDNC (ASDIS,
DSCRD, ASSMB, TSI_G), SPAN, DISC_MAPQ, NUMPARTS (assembly parts), NM, MAPQ.
FORMAT per sample: SR, DR, AD, DP, GQ.
"""
from __future__ import annotations

from ..model import BreakpointCall, Evidence
from .base import iter_records, parse_info, parse_format_sample, parse_bnd_alt, normalise_order


def parse_svaba(path: str, sample: str) -> list[BreakpointCall]:
    seen_mates: set[str] = set()
    out: list[BreakpointCall] = []
    for f in iter_records(path):
        chrom, pos, vid, ref, alt, qual, filt, info_str = f[:8]
        info = parse_info(info_str)
        if "SVTYPE" not in info:
            continue
        sv_type = info.get("SVTYPE", "BND")
        mate_id = info.get("MATEID", "")
        if vid in seen_mates or mate_id in seen_mates:
            continue
        seen_mates.add(vid)
        precise = "IMPRECISE" not in info
        pass_filter = (filt == "PASS")
        try:
            q = float(qual) if qual not in ("", ".") else 0.0
        except ValueError:
            q = 0.0

        parsed = parse_bnd_alt(alt)
        if parsed:
            mchrom, mpos, sa, sb = parsed
            chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
                chrom, int(pos), sa, mchrom, mpos, sb
            )
        else:
            end = int(info.get("END", pos))
            chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
                chrom, int(pos), "+", chrom, end, "-"
            )

        evdnc = info.get("EVDNC", "")  # ASDIS | DSCRD | ASSMB | TSI_G
        numparts = int(info.get("NUMPARTS") or 0)
        mq = int(info.get("MAPQ") or info.get("DISC_MAPQ") or 0)
        fmt_sample = {}
        if len(f) >= 10:
            fmt_sample = parse_format_sample(f[8], f[9])
        sr = int(_first_numeric(fmt_sample.get("SR")))
        dr = int(_first_numeric(fmt_sample.get("DR")))
        ad = int(_first_numeric(fmt_sample.get("AD")))
        dp = int(_first_numeric(fmt_sample.get("DP")))
        vaf = (ad / dp) if dp > 0 else 0.0

        # If EVDNC indicates DSCRD only, treat sr=0; if ASSMB only, treat as assembly+split combined
        if evdnc == "DSCRD":
            sr_evi, dr_evi = 0, max(dr, 1)
        elif evdnc == "ASSMB":
            sr_evi, dr_evi = max(sr, 1), 0
        else:  # ASDIS or unknown — both
            sr_evi, dr_evi = sr, dr

        ev = Evidence(
            caller="svaba",
            split_reads=sr_evi,
            discordant_pairs=dr_evi,
            assembly_contigs=numparts if numparts > 0 else (1 if evdnc in ("ASDIS", "ASSMB") else 0),
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
            record_id=vid,
        ))
    return out


def _first_numeric(s: str | None) -> int:
    if s is None or s in ("", "."):
        return 0
    head = s.split(",")[0]
    try:
        return int(head)
    except ValueError:
        try:
            return int(float(head))
        except ValueError:
            return 0
