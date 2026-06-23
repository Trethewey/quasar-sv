#!/usr/bin/env python3
"""For each IG / TR locus in a sample, count the distribution of mate
chromosomes and (when on a driver chromosome) mate-locus annotation.

This directly answers: at the IGH/IGL/IGK loci, what is the dominant
non-self mate chromosome? If IGH ↔ chr3:BCL6 is real, IGH should have a
big chr3 bump that IGL and IGK do not.
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
    ("CIITA",   "16",  10_971_055,  11_018_791),
    ("ALK",      "2",  29_192_774,  29_921_566),
    ("CCND3",    "6",  41_934_935,  41_956_141),
    ("IRF4",     "6",     391_739,     411_447),
    ("BCL3",    "19",  44_747_756,  44_760_044),
    ("BCL11A",   "2",  60_450_519,  60_553_544),
    ("EZH2",     "7", 148_807_374, 148_884_326),
    ("MEF2B",   "19",  19_143_719,  19_170_610),
    ("BCL10",    "1",  85_458_264,  85_466_361),
]


def _strip(c: str) -> str:
    return c[3:] if c.startswith("chr") else c


def _locus_for(chrom: str, pos: int, table) -> str:
    chrom = _strip(chrom)
    for name, c, s, e in table:
        if chrom == c and s <= pos <= e:
            return name
    return ""


def scan(cram, reference, sample, max_reads_per_locus=2_000_000,
         min_mapq=20, drop_artefact=True):
    """For each IG locus: count mate chrom distribution + per-driver mate counts.

    Returns dict[locus_name][category] = Counter.
    """
    open_kwargs = {"reference_filename": reference} if cram.endswith(".cram") else {}
    mode = "rc" if cram.endswith(".cram") else "rb"
    sam = pysam.AlignmentFile(cram, mode, **open_kwargs)
    refs = set(sam.references)
    chr_prefix = any(r.startswith("chr") for r in refs)
    def to_ref(c):
        return (c if c.startswith("chr") else f"chr{c}") if chr_prefix else _strip(c)

    art_chrom, art_start, art_end = "2", 32_915_800, 32_916_800
    out: dict[str, dict] = {}
    for ig_name, ig_chrom, ig_start, ig_end in IG_LOCI:
        ig_ref = to_ref(ig_chrom)
        if ig_ref not in refs:
            continue
        mate_chrom = Counter()
        driver_mates = Counter()
        sa_chrom = Counter()
        sa_driver = Counter()
        n_reads = 0
        for read in sam.fetch(ig_ref, ig_start, ig_end):
            n_reads += 1
            if n_reads > max_reads_per_locus:
                break
            if read.is_unmapped or read.is_secondary or read.is_duplicate:
                continue
            if read.mapping_quality < min_mapq:
                continue
            # Mate-based PE evidence
            if read.is_paired and not read.mate_is_unmapped:
                m_chrom = sam.get_reference_name(read.next_reference_id) or ""
                m_pos = read.next_reference_start
                m_chrom_n = _strip(m_chrom)
                # Discard intra-locus
                if m_chrom_n == ig_chrom and ig_start - 5_000 <= m_pos <= ig_end + 5_000:
                    pass
                else:
                    if drop_artefact and m_chrom_n == art_chrom and art_start <= m_pos <= art_end:
                        pass
                    else:
                        mate_chrom[m_chrom_n] += 1
                        drv = _locus_for(m_chrom, m_pos, DRIVERS)
                        if drv:
                            driver_mates[drv] += 1
            # SA-tag based SR evidence
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
                    s_chrom_n = _strip(s_chrom)
                    if s_chrom_n == ig_chrom and ig_start - 5_000 <= s_pos <= ig_end + 5_000:
                        pass
                    elif drop_artefact and s_chrom_n == art_chrom and art_start <= s_pos <= art_end:
                        pass
                    else:
                        sa_chrom[s_chrom_n] += 1
                        drv = _locus_for(s_chrom, s_pos, DRIVERS)
                        if drv:
                            sa_driver[drv] += 1
        out[ig_name] = dict(
            mate_chrom=mate_chrom,
            driver_mates=driver_mates,
            sa_chrom=sa_chrom,
            sa_driver=sa_driver,
            n_reads=n_reads,
        )
    sam.close()
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reference", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("samples", nargs="+", help="sample=cram pairs")
    args = p.parse_args()

    rows = []
    for spec in args.samples:
        name, cram = spec.split("=", 1)
        print(f"[scan] {name}: {cram}", file=sys.stderr)
        d = scan(cram, args.reference, name)
        print(f"\n=== {name} ===")
        for ig in sorted(d):
            r = d[ig]
            print(f"  {ig}: n_reads={r['n_reads']}")
            top_chr = ", ".join(f"chr{c}={n}" for c, n in r["mate_chrom"].most_common(5))
            print(f"    mate chrom top 5:    {top_chr}")
            top_drv = ", ".join(f"{c}={n}" for c, n in r["driver_mates"].most_common(5))
            print(f"    driver mate hits:    {top_drv if top_drv else '(none)'}")
            top_sa_chr = ", ".join(f"chr{c}={n}" for c, n in r["sa_chrom"].most_common(5))
            print(f"    SA chrom top 5:      {top_sa_chr}")
            top_sa_drv = ", ".join(f"{c}={n}" for c, n in r["sa_driver"].most_common(5))
            print(f"    driver SA hits:      {top_sa_drv if top_sa_drv else '(none)'}")
        for ig, r in d.items():
            for c, n in r["mate_chrom"].most_common(10):
                rows.append((name, ig, "mate_chrom", c, n))
            for drv, n in r["driver_mates"].most_common(20):
                rows.append((name, ig, "driver_mate", drv, n))
            for c, n in r["sa_chrom"].most_common(10):
                rows.append((name, ig, "sa_chrom", c, n))
            for drv, n in r["sa_driver"].most_common(20):
                rows.append((name, ig, "sa_driver", drv, n))

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "ig_locus", "category", "key", "count"])
        for row in rows:
            w.writerow(row)


if __name__ == "__main__":
    main()
