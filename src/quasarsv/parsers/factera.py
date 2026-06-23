"""FACTERA parser.

FACTERA is panel-only and writes its own tabular `.fusions.txt` output.
Schema (tab-separated): Est_Type Region1 Region2 Break1 Break2 Break_support1
Break_support2 Break_offset Orientation Order Both_exons Non-templated_seq
size_(kb) Coding_potential Total_depth ...

Treat as: split-read support per breakpoint = Break_support1 / Break_support2.
Two breakpoints define the fusion; positions in Break1/Break2.
"""
from __future__ import annotations

from pathlib import Path

from ..model import BreakpointCall, Evidence
from .base import normalise_order


def parse_factera(path: str, sample: str) -> list[BreakpointCall]:
    p = Path(path)
    if not p.exists():
        return []
    out: list[BreakpointCall] = []
    with open(p, encoding="utf-8") as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            cols = line.split("\t")
            if header is None:
                header = [c.strip() for c in cols]
                continue
            row = dict(zip(header, cols))
            try:
                # Break1/Break2 are like 'chr6:117658326'; Region1/Region2 are gene symbols
                ca, pa = row["Break1"].split(":")
                cb, pb = row["Break2"].split(":")
                chrom_a = _strip_chr(ca)
                chrom_b = _strip_chr(cb)
                pos_a = int(pa)
                pos_b = int(pb)
            except (KeyError, ValueError):
                continue
            sr_a = _safe_int(row.get("Break_support1"))
            sr_b = _safe_int(row.get("Break_support2"))
            sr = min(sr_a, sr_b) if (sr_a and sr_b) else max(sr_a, sr_b)
            orient = row.get("Orientation", "")
            sa, sb = _orient_to_strands(orient)
            chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b = normalise_order(
                chrom_a, pos_a, sa, chrom_b, pos_b, sb
            )
            ev = Evidence(
                caller="factera",
                split_reads=sr,
                discordant_pairs=0,
                assembly_contigs=1,            # FACTERA contigs a breakpoint by design
                mapq=0,
                vaf=0.0,
                filter_pass=True,              # FACTERA only writes passing fusions
                precise=True,
                raw_qual=float(sr),
            )
            out.append(BreakpointCall(
                sample=sample, chrom_a=chrom_a, pos_a=pos_a, strand_a=strand_a,
                chrom_b=chrom_b, pos_b=pos_b, strand_b=strand_b,
                sv_type="BND", evidence=ev,
                record_id=f"factera_{chrom_a}_{pos_a}_{chrom_b}_{pos_b}",
            ))
    return out


def _strip_chr(s: str) -> str:
    return s[3:] if s.lower().startswith("chr") else s


def _safe_int(s: str | None) -> int:
    if s is None or s == "":
        return 0
    try:
        return int(s)
    except ValueError:
        return 0


def _orient_to_strands(orient: str) -> tuple[str, str]:
    o = orient.strip()
    if o in ("HT", "+-"):
        return "+", "-"
    if o in ("TH", "-+"):
        return "-", "+"
    if o in ("HH", "++"):
        return "+", "+"
    if o in ("TT", "--"):
        return "-", "-"
    return "+", "-"
