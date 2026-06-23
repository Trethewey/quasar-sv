"""Shared VCF helpers — open (gzip-aware), parse BND ALT, normalise breakpoint order."""
from __future__ import annotations

import gzip
import re
from typing import Iterator

BND_ALT_RE = re.compile(r"([ACGTN\.]*)([\[\]])([^:\[\]]+):(\d+)([\[\]])([ACGTN\.]*)", re.IGNORECASE)


def open_vcf(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def iter_records(path: str) -> Iterator[list[str]]:
    """Yield non-header VCF records split on tab."""
    with open_vcf(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            yield line.split("\t")


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in info.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = "1"  # flag
    return out


def parse_format_sample(format_field: str, sample_field: str) -> dict[str, str]:
    keys = format_field.split(":")
    vals = sample_field.split(":")
    return dict(zip(keys, vals))


def parse_bnd_alt(alt: str) -> tuple[str, int, str, str] | None:
    """Parse VCF BND ALT (e.g. 'N[chr2:32916489[') -> (mate_chrom, mate_pos, strand_a, strand_b).

    Strand convention per VCF 4.3:
      t[p[  joined after t  -> strand_a='+', mate left of break  -> strand_b='+'
      t]p]  joined after t  -> strand_a='+', mate right of break -> strand_b='-'
      ]p]t  joined before t -> strand_a='-', mate right of break -> strand_b='-'
      [p[t  joined before t -> strand_a='-', mate left of break  -> strand_b='+'
    """
    m = BND_ALT_RE.search(alt)
    if not m:
        return None
    pre, b1, mchrom, mpos, b2, post = m.groups()
    if b1 != b2:
        return None
    # bracket orientation determines strands
    bracket = b1
    if pre and not post:           # t[p[ or t]p]
        strand_a = "+"
        strand_b = "+" if bracket == "[" else "-"
    elif post and not pre:         # [p[t or ]p]t
        strand_a = "-"
        strand_b = "+" if bracket == "[" else "-"
    else:
        return None
    return mchrom, int(mpos), strand_a, strand_b


def normalise_order(chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b):
    """Return breakpoints in a deterministic order so two callers describing the
    same translocation produce the same (chrom_a, pos_a)."""
    a = (_chrom_key(chrom_a), pos_a)
    b = (_chrom_key(chrom_b), pos_b)
    if a <= b:
        return chrom_a, pos_a, strand_a, chrom_b, pos_b, strand_b
    return chrom_b, pos_b, strand_b, chrom_a, pos_a, strand_a


def _chrom_key(c: str) -> tuple[int, str]:
    s = c[3:] if c.lower().startswith("chr") else c
    try:
        return (0, f"{int(s):02d}")
    except ValueError:
        return (1, s)
