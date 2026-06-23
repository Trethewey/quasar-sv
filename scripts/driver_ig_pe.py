#!/usr/bin/env python3
"""Bidirectional driver↔IG PE/SR support probe.

For each driver locus, scan reads and count mate / SA distributions, with
special attention to mates that land in each candidate IG locus (and the
chr2:32916 polyG artefact, which absorbs mis-routed IG-switch reads).

This gives a direct discriminator for "which IG is the partner of this driver"
that does not depend on the polyG-rescue heuristic.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import pysam

IG_LOCI = [
    ("IGH",       "14", 105_583_730, 106_879_812),
    ("IGK",       "2",   88_857_361,  90_235_368),
    ("IGL",      "22",   22_026_076,  22_922_913),
    ("TRA",      "14",   21_621_904,  22_552_132),
    ("TRB",       "7",  142_299_011, 142_813_287),
    ("TRG",       "7",   38_240_024,  38_368_055),
    ("TRD",      "14",   22_422_546,  22_466_577),
]
DRIVERS = [
    ("BCL6",     "3", 187_420_000, 187_800_000),
    ("BCL2",    "18",  63_123_346,  63_320_128),
    ("MYC",      "8", 127_735_434, 127_741_434),
    ("CCND1",   "11",  69_641_087,  69_654_474),
    ("MALT1",   "18",  58_671_498,  58_754_525),
    ("PAX5",     "9",  36_833_275,  37_034_103),
    ("FOXP1",    "3",  71_247_033,  71_633_153),
    ("BCL3",    "19",  44_747_756,  44_760_044),
    ("BCL11A",   "2",  60_450_519,  60_553_544),
    ("IRF4",     "6",     391_739,     411_447),
    ("EZH2",     "7", 148_807_374, 148_884_326),
]
ART = ("2", 32_915_800, 32_916_800)


def _strip(c):
    return c[3:] if c.startswith("chr") else c


def _in(table, chrom, pos):
    c = _strip(chrom)
    for name, ch, s, e in table:
        if c == ch and s <= pos <= e:
            return name
    return ""


def scan_driver(sam, chr_prefix, drv_name, drv_chrom, drv_start, drv_end,
                min_mapq=20, max_reads=2_000_000):
    ref = (drv_chrom if drv_chrom.startswith("chr") else f"chr{drv_chrom}") if chr_prefix else _strip(drv_chrom)
    if ref not in sam.references:
        return None
    mate_chrom = Counter()
    mate_ig = Counter()      # mate landing in any IG locus
    sa_chrom = Counter()
    sa_ig = Counter()
    mate_artefact = 0
    sa_artefact = 0
    n_reads = 0
    n_chim = 0
    for read in sam.fetch(ref, drv_start, drv_end):
        n_reads += 1
        if n_reads > max_reads:
            break
        if read.is_unmapped or read.is_secondary or read.is_duplicate:
            continue
        if read.mapping_quality < min_mapq:
            continue
        # PE: mate elsewhere
        if read.is_paired and not read.mate_is_unmapped:
            m_chrom = sam.get_reference_name(read.next_reference_id) or ""
            m_pos = read.next_reference_start
            mn = _strip(m_chrom)
            same_locus = (mn == drv_chrom and drv_start - 5_000 <= m_pos <= drv_end + 5_000)
            if not same_locus:
                n_chim += 1
                ig_hit = _in(IG_LOCI, m_chrom, m_pos)
                if ig_hit:
                    mate_ig[ig_hit] += 1
                if mn == ART[0] and ART[1] <= m_pos <= ART[2]:
                    mate_artefact += 1
                else:
                    mate_chrom[mn] += 1
        # SR: SA tag
        try:
            sa = read.get_tag("SA")
        except KeyError:
            sa = ""
        if sa:
            first = sa.split(";", 1)[0]
            parts = first.split(",")
            if len(parts) >= 2:
                s_chrom = parts[0]
                s_pos = int(parts[1])
                sn = _strip(s_chrom)
                same_locus = (sn == drv_chrom and drv_start - 5_000 <= s_pos <= drv_end + 5_000)
                if not same_locus:
                    ig_hit = _in(IG_LOCI, s_chrom, s_pos)
                    if ig_hit:
                        sa_ig[ig_hit] += 1
                    if sn == ART[0] and ART[1] <= s_pos <= ART[2]:
                        sa_artefact += 1
                    else:
                        sa_chrom[sn] += 1
    return dict(
        driver=drv_name,
        n_reads=n_reads,
        n_chim=n_chim,
        mate_chrom=mate_chrom,
        mate_ig=mate_ig,
        mate_artefact=mate_artefact,
        sa_chrom=sa_chrom,
        sa_ig=sa_ig,
        sa_artefact=sa_artefact,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reference", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("samples", nargs="+", help="sample=cram pairs")
    args = p.parse_args()

    rows = []
    for spec in args.samples:
        name, cram = spec.split("=", 1)
        print(f"[scan] {name}", file=sys.stderr)
        open_kwargs = {"reference_filename": args.reference} if cram.endswith(".cram") else {}
        mode = "rc" if cram.endswith(".cram") else "rb"
        sam = pysam.AlignmentFile(cram, mode, **open_kwargs)
        chr_prefix = any(r.startswith("chr") for r in sam.references)
        print(f"\n=== {name} ===")
        for drv_name, drv_chrom, drv_start, drv_end in DRIVERS:
            r = scan_driver(sam, chr_prefix, drv_name, drv_chrom, drv_start, drv_end)
            if r is None:
                continue
            top_ig_mate = sorted(r["mate_ig"].items(), key=lambda kv: -kv[1])
            top_ig_sa = sorted(r["sa_ig"].items(), key=lambda kv: -kv[1])
            top_mate = ", ".join(f"chr{c}={n}" for c, n in r["mate_chrom"].most_common(5))
            print(f"  {drv_name:<8} reads={r['n_reads']:>8}  chimeric={r['n_chim']:>8}  "
                  f"mate@artefact={r['mate_artefact']:>5}  sa@artefact={r['sa_artefact']:>5}")
            print(f"     mate-in-IG: {top_ig_mate or '(none)'}")
            print(f"     sa-in-IG:   {top_ig_sa or '(none)'}")
            print(f"     mate top chr (excl. art): {top_mate}")
            for ig, n in r["mate_ig"].items():
                rows.append((name, drv_name, "mate_ig", ig, n))
            for ig, n in r["sa_ig"].items():
                rows.append((name, drv_name, "sa_ig", ig, n))
            rows.append((name, drv_name, "mate_artefact", "chr2_32916", r["mate_artefact"]))
            rows.append((name, drv_name, "sa_artefact", "chr2_32916", r["sa_artefact"]))
        sam.close()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "driver", "category", "key", "count"])
        for row in rows:
            w.writerow(row)


if __name__ == "__main__":
    main()
