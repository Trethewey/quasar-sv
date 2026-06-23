#!/usr/bin/env python3
"""Emit a BED of lymphoma driver + IG/TR + artefact loci, padded for caller use.

Manta/GRIDSS/Delly/SvABA/TIDDIT all accept a target-region BED to restrict the
call space — running on the full WGS would take 6-30 hours per CRAM but
targeted to ~50 Mb of loci it completes in 10-30 minutes per tool per CRAM.

Default pad of 200 kb each side so we don't lose breakpoints near the locus
boundary or chimeric mates landing just outside it.
"""
from __future__ import annotations

import argparse
import csv
import sys
from importlib.resources import files as pkg_files
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pad-bp", type=int, default=200_000,
                   help="bp padding either side of each locus (default 200,000)")
    p.add_argument("--chr-prefix", action="store_true", default=True,
                   help="emit chr-prefixed contigs to match GRCh38 with full decoys")
    p.add_argument("--output", required=True, help="output BED path")
    args = p.parse_args()

    sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
    from quasarsv.annotate import load_builtin_loci
    from quasarsv.qc import _load_artefact_loci

    rows: list[tuple[str, int, int, str]] = []
    for g in load_builtin_loci():
        chrom = f"chr{g.chrom}" if args.chr_prefix else g.chrom
        start = max(0, g.start - args.pad_bp)
        end = g.end + args.pad_bp
        rows.append((chrom, start, end, f"{g.gene}|{g.role}"))
    for chrom_n, start, end, notes in _load_artefact_loci():
        chrom = f"chr{chrom_n}" if args.chr_prefix else chrom_n
        rows.append((chrom, max(0, start - args.pad_bp), end + args.pad_bp,
                     f"artefact|{notes[:40]}"))

    # Sort + merge overlaps
    rows.sort(key=lambda r: (r[0], r[1]))
    merged: list[tuple[str, int, int, str]] = []
    for chrom, start, end, name in rows:
        if merged and merged[-1][0] == chrom and merged[-1][2] >= start:
            prev = merged[-1]
            merged[-1] = (chrom, prev[1], max(prev[2], end), prev[3] + "," + name)
        else:
            merged.append((chrom, start, end, name))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for r in merged:
            w.writerow(r)
    total_bp = sum(e - s for _, s, e, _ in merged)
    print(f"[bed] {len(merged)} regions, {total_bp/1e6:.1f} Mb total -> {out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
