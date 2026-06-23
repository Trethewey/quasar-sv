#!/usr/bin/env python3
"""SA-tag chromosome distribution at the chr2:32916 polyG artefact.

For one or more CRAMs, scans the artefact locus and tallies SA-tag mate
chromosomes. Also tallies per-chromosome how many SA tags land inside any
known IG / driver locus window (so we can see, for example, whether chr14
SAs are dominantly inside the IGH window or the TRA window).
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

ART_CHROM = "chr2"
ART_START = 32_915_800
ART_END = 32_916_800

# IG / TR locus windows on each lineage's chromosome
IG_LOCI = [
    ("IGH",       "14", 105_583_730, 106_879_812),
    ("IGK",       "2",   88_857_361,  90_235_368),
    ("IGL",      "22",   22_026_076,  22_922_913),
    ("TRA",      "14",   21_621_904,  22_552_132),
    ("TRB",       "7",  142_299_011, 142_813_287),
    ("TRG",       "7",   38_240_024,  38_368_055),
    ("TRD",      "14",   22_422_546,  22_466_577),
]
DRIVERS_OF_INTEREST = [
    ("BCL6",     "3", 187_420_000, 187_800_000),
    ("BCL2",    "18",  63_123_346,  63_320_128),
    ("MYC",      "8", 127_735_434, 127_741_434),
    ("CCND1",   "11",  69_641_087,  69_654_474),
    ("MALT1",   "18",  58_671_498,  58_754_525),
]


def _strip(c: str) -> str:
    return c[3:] if c.startswith("chr") else c


def _in_window(chrom: str, pos: int, windows: list[tuple[str, str, int, int]]) -> str:
    chrom = _strip(chrom)
    for name, c, s, e in windows:
        if chrom == c and s <= pos <= e:
            return name
    return ""


def scan(cram: str, reference: str, sample: str,
         max_reads: int = 5_000_000, apply_qc_filters: bool = True,
         use_mate: bool = False) -> dict:
    open_kwargs = {"reference_filename": reference} if cram.endswith(".cram") else {}
    mode = "rc" if cram.endswith(".cram") else "rb"
    sam = pysam.AlignmentFile(cram, mode, **open_kwargs)
    refs = set(sam.references)
    chr_prefix = any(r.startswith("chr") for r in refs)
    art = ART_CHROM if chr_prefix else _strip(ART_CHROM)

    chrom_counts = Counter()
    locus_counts = Counter()
    chr14_position_buckets = Counter()
    chr22_position_buckets = Counter()
    chr2_position_buckets = Counter()
    n_reads = 0
    n_with_sa = 0
    n_with_distant_mate = 0
    for read in sam.fetch(art, ART_START, ART_END):
        n_reads += 1
        if n_reads > max_reads:
            break
        if read.is_unmapped:
            continue
        if apply_qc_filters:
            if read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality < 20:
                continue
        try:
            sa = read.get_tag("SA")
            has_sa = True
        except KeyError:
            has_sa = False
        if has_sa:
            n_with_sa += 1
            first = sa.split(";", 1)[0]
            parts = first.split(",")
            if len(parts) >= 2:
                sa_chrom = parts[0]
                sa_pos = int(parts[1])
                chrom_counts[_strip(sa_chrom)] += 1
                loc_name = _in_window(sa_chrom, sa_pos, IG_LOCI + DRIVERS_OF_INTEREST)
                if loc_name:
                    locus_counts[loc_name] += 1
                sa_n = _strip(sa_chrom)
                bucket = sa_pos // 1_000_000
                if sa_n == "14":
                    chr14_position_buckets[bucket] += 1
                elif sa_n == "22":
                    chr22_position_buckets[bucket] += 1
                elif sa_n == "2":
                    chr2_position_buckets[bucket] += 1
        if use_mate and read.is_paired and not read.mate_is_unmapped:
            mate_chrom = sam.get_reference_name(read.next_reference_id) or ""
            mate_pos = read.next_reference_start
            mate_n = _strip(mate_chrom)
            # Skip intra-artefact PE
            if mate_n == _strip(ART_CHROM) and ART_START - 1000 <= mate_pos <= ART_END + 1000:
                continue
            n_with_distant_mate += 1
            chrom_counts[mate_n + "_mate"] += 1
            loc_name = _in_window(mate_chrom, mate_pos, IG_LOCI + DRIVERS_OF_INTEREST)
            if loc_name:
                locus_counts[loc_name + "_mate"] += 1
            bucket = mate_pos // 1_000_000
            if mate_n == "14":
                chr14_position_buckets[("mate", bucket)] += 1
            elif mate_n == "22":
                chr22_position_buckets[("mate", bucket)] += 1
    sam.close()
    chrom_counts["__n_with_distant_mate"] = n_with_distant_mate
    return {
        "sample": sample,
        "n_reads": n_reads,
        "n_with_sa": n_with_sa,
        "chrom_counts": chrom_counts.most_common(),
        "locus_counts": locus_counts.most_common(),
        "chr14_position_buckets": chr14_position_buckets.most_common(8),
        "chr22_position_buckets": chr22_position_buckets.most_common(8),
        "chr2_position_buckets": chr2_position_buckets.most_common(8),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reference", required=True)
    p.add_argument("--out", required=True,
                   help="TSV: sample, chrom or locus, count")
    p.add_argument("--no-filters", action="store_true",
                   help="skip mapq/duplicate/secondary filters (match handover bash diagnostic)")
    p.add_argument("--with-mate", action="store_true",
                   help="also tally MATE chromosome / locus distributions")
    p.add_argument("--max-reads", type=int, default=5_000_000)
    p.add_argument("samples", nargs="+",
                   help="sample=cram pairs, e.g. Karpas1106P=/path/to.cram")
    args = p.parse_args()

    rows = []
    summaries = []
    for spec in args.samples:
        if "=" not in spec:
            sys.exit(f"sample spec must be 'name=cram_path': {spec}")
        name, cram = spec.split("=", 1)
        print(f"[scan] {name}: {cram}", file=sys.stderr)
        r = scan(cram, args.reference, name, max_reads=args.max_reads,
                 apply_qc_filters=not args.no_filters, use_mate=args.with_mate)
        summaries.append(r)
        for chrom, n in r["chrom_counts"]:
            rows.append((name, "chrom", chrom, n))
        for loc, n in r["locus_counts"]:
            rows.append((name, "locus", loc, n))
        for b, n in r["chr14_position_buckets"]:
            rows.append((name, "chr14_mb_bucket", str(b), n))
        for b, n in r["chr22_position_buckets"]:
            rows.append((name, "chr22_mb_bucket", str(b), n))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "category", "key", "count"])
        for row in rows:
            w.writerow(row)

    for r in summaries:
        print(f"\n=== {r['sample']} (n_reads={r['n_reads']}, with_sa={r['n_with_sa']}) ===")
        print("Top chrom SA targets:")
        for c, n in r["chrom_counts"][:15]:
            print(f"  {c:<14} {n:>8}")
        print("Top locus SA / mate targets:")
        for loc, n in r["locus_counts"][:15]:
            print(f"  {loc:<14} {n:>8}")
        print("Top chr14 Mb buckets:")
        for b, n in r["chr14_position_buckets"]:
            print(f"  chr14:{b} {n:>8}")
        print("Top chr22 Mb buckets:")
        for b, n in r["chr22_position_buckets"]:
            print(f"  chr22:{b} {n:>8}")


if __name__ == "__main__":
    main()
