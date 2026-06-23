"""GRIDSS2 VCF parser.

GRIDSS emits BND-only records, every breakpoint paired via MATEID (EVENT
identifies the SV; one record per end). Per-end evidence in FORMAT and INFO:
  SR/SRQ            split reads (single break)
  RP/RPQ            discordant read pairs
  ASRP/ASSR         pairs/reads contributing through assembly
  IC/ASC            assembly contig count
  BSC/BSCQ          soft-clip count
  AF                allele fraction
"""
from __future__ import annotations

from ..model import BreakpointCall, Evidence
from .base import iter_records, parse_info, parse_format_sample, parse_bnd_alt, normalise_order


def parse_gridss(path: str, sample: str) -> list[BreakpointCall]:
    seen_events: set[str] = set()
    out: list[BreakpointCall] = []

    for f in iter_records(path):
        chrom, pos, vid, ref, alt, qual, filt, info_str = f[:8]
        info = parse_info(info_str)
        if info.get("SVTYPE") != "BND":
            continue
        event = info.get("EVENT", vid)
        if event in seen_events:
            continue
        seen_events.add(event)

        parsed = parse_bnd_alt(alt)
        if parsed is None:
            continue
        mchrom, mpos, sa, sb = parsed
        chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
            chrom, int(pos), sa, mchrom, mpos, sb
        )

        fmt_sample = {}
        if len(f) >= 10:
            fmt_sample = parse_format_sample(f[8], f[9])

        sr = _i(fmt_sample.get("SR") or info.get("SR"))
        sr_total = sr + _i(fmt_sample.get("ASSR") or info.get("ASSR"))
        rp = _i(fmt_sample.get("RP") or info.get("RP"))
        rp_total = rp + _i(fmt_sample.get("ASRP") or info.get("ASRP"))
        ic = _i(info.get("IC"))
        asm_contigs = ic + (1 if _i(info.get("AS")) > 0 else 0)
        bsc = _i(info.get("BSC"))
        af = _f(fmt_sample.get("AF") or info.get("AF"))
        mq = _i(info.get("MQ"))
        try:
            q = float(qual) if qual not in ("", ".") else 0.0
        except ValueError:
            q = 0.0
        pass_filter = (filt == "PASS")
        precise = "IMPRECISE" not in info

        ci = _parse_ci(info.get("CIPOS"))
        cir = _parse_ci(info.get("CIRPOS") or info.get("CIPOS"))

        ev = Evidence(
            caller="gridss",
            split_reads=sr_total,
            discordant_pairs=rp_total,
            assembly_contigs=asm_contigs,
            soft_clips=bsc,
            mapq=mq,
            vaf=af,
            filter_pass=pass_filter,
            precise=precise,
            raw_qual=q,
        )
        out.append(BreakpointCall(
            sample=sample, chrom_a=chrom_a, pos_a=pos_a, strand_a=strand_a,
            chrom_b=chrom_b, pos_b=pos_b, strand_b=strand_b,
            sv_type="BND", evidence=ev,
            cipos_a=ci, cipos_b=cir, record_id=event,
        ))
    return out


def _i(s) -> int:
    if s is None or s == "" or s == ".":
        return 0
    try:
        return int(s)
    except ValueError:
        try:
            return int(float(s))
        except ValueError:
            return 0


def _f(s) -> float:
    if s is None or s == "" or s == ".":
        return 0.0
    try:
        return float(s)
    except ValueError:
        return 0.0


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
