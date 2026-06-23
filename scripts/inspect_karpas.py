"""Inspect Karpas-1106P scan output for the expected BCL6-IGH t(3;14)."""
import sys
sys.path.insert(0, "src")
from quasarsv.model import read_fusion_calls_tsv

calls = read_fusion_calls_tsv("output/wgs_karpas1106p_scan/ERR9188549_Karpas1106P.fusions.tsv")
print(f"total candidates: {len(calls)}")

def in_locus(c, side, chrom, lo, hi):
    return getattr(c, f"chrom_{side}") == chrom and lo <= getattr(c, f"pos_{side}") <= hi

# BCL6 (chr3:187.7M) <-> IGH (chr14:105.5-107M)
bcl6 = lambda c, s: in_locus(c, s, "3", 187_700_000, 187_750_000)
igh = lambda c, s: in_locus(c, s, "14", 105_500_000, 107_000_000)

bcl6_igh = [c for c in calls if (bcl6(c, "a") and igh(c, "b")) or (bcl6(c, "b") and igh(c, "a"))]
print(f"\nBCL6 <-> IGH (t(3;14)) candidates: {len(bcl6_igh)}")
bcl6_igh.sort(key=lambda c: -(c.split_reads + c.discordant_pairs))
for c in bcl6_igh[:8]:
    print(f"  {c.chrom_a}:{c.pos_a:>10}{c.strand_a} <-> {c.chrom_b}:{c.pos_b:>10}{c.strand_b}  "
          f"gA={c.gene_a or '-':<8} gB={c.gene_b or '-':<8} SR={c.split_reads:>4} PE={c.discordant_pairs:>4} "
          f"tier={c.tier} known={c.known_partner}")

# Any call touching IGH
igh_calls = [c for c in calls if igh(c, "a") or igh(c, "b")]
print(f"\nany call touching IGH (chr14:105.5-107M): {len(igh_calls)}")
igh_calls.sort(key=lambda c: -(c.split_reads + c.discordant_pairs))
for c in igh_calls[:10]:
    print(f"  {c.chrom_a}:{c.pos_a:>10}{c.strand_a} <-> {c.chrom_b}:{c.pos_b:>10}{c.strand_b}  "
          f"gA={c.gene_a or '-':<8} gB={c.gene_b or '-':<8} driver={c.driver_locus:<12} "
          f"SR={c.split_reads:>4} PE={c.discordant_pairs:>4}")

# Inter-chr calls with highest support overall
inter = [c for c in calls if c.chrom_a != c.chrom_b]
inter.sort(key=lambda c: -(c.split_reads + c.discordant_pairs))
print(f"\nTop 10 inter-chr calls by support:")
for c in inter[:10]:
    print(f"  {c.chrom_a}:{c.pos_a:>10} <-> {c.chrom_b}:{c.pos_b:>10}  "
          f"gA={c.gene_a or '-':<8} gB={c.gene_b or '-':<8} SR={c.split_reads:>4} PE={c.discordant_pairs:>4}")
