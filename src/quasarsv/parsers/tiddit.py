"""TIDDIT VCF parser.

TIDDIT 3.x emits INV, DEL, DUP, INS, TDUP, BND. Evidence lives in two places and
the distinction matters — reading it wrongly silently zeroes every call:

  FORMAT (per-sample, authoritative):
    DV  discordant pairs supporting the variant
    RV  split ("reference-variant") reads supporting the variant
  INFO (Number=2, "read-pairs,split-reads" — NOT scalars):
    LTE  links to the event         e.g. LTE=26,5
    LFA  links from window A        e.g. LFA=30,6
    LFB  links from window B

An earlier version read ``DV``/``RV`` out of INFO (they are FORMAT), passed the
comma-pair ``LFA``/``LFB`` through an int() that threw and defaulted to 0, and
summed ``LTA``/``LTB``, which TIDDIT does not emit at all (the event field is
``LTE``). Every record therefore scored split_reads=0 / discordant_pairs=0, fell
below the single-caller tier thresholds, and landed at T3 — producing a
benchmark result of tp=0 AND fp=0, i.e. TIDDIT appearing to detect nothing
whatsoever. That was our bug being reported as TIDDIT's performance.

Quality is in QUAL (per-record). FILTER is PASS / various tiddit-specific tags.
"""
from __future__ import annotations

from ..model import BreakpointCall, Evidence
from .base import (
    iter_records, parse_info, parse_bnd_alt, normalise_order,
    parse_format_sample,
)


def _int(v, default=0):
    if v is None or v == "" or v == ".":
        return default
    try:
        return int(v)
    except ValueError:
        return default


def _pair(v, idx: int, default=0):
    """Read one element of a TIDDIT Number=2 'read-pairs,split-reads' field."""
    if not v or v == ".":
        return default
    parts = str(v).split(",")
    if idx >= len(parts):
        return default
    return _int(parts[idx], default)


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

        # Prefer the per-sample FORMAT counts; fall back to INFO/LTE (the
        # "links to event" pair) when a VCF carries no sample column.
        fmt_sample = {}
        if len(f) >= 10:
            fmt_sample = parse_format_sample(f[8], f[9])
        disc = _int(fmt_sample.get("DV"), _pair(info.get("LTE"), 0))
        split = _int(fmt_sample.get("RV"), _pair(info.get("LTE"), 1))

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
