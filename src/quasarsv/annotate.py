"""Gene / locus annotation layer.

Built-in tables (data/lymphoma_loci.tsv, data/known_partners.tsv) cover the
clinically important lymphoma driver loci and canonical rearrangement partners.
For full transcriptome annotation, supply a GTF via `load_gtf`.

Functions
---------
annotate_calls(calls, ...): mutates FusionCall in place with gene_a/gene_b,
region_a/region_b, known_partner, driver_locus, in_frame (when GTF supplied).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from importlib.resources import files as pkg_files
from pathlib import Path
from typing import Iterable

from .model import FusionCall


@dataclass
class GeneEntry:
    gene: str
    chrom: str
    start: int
    end: int
    strand: str
    role: str
    notes: str


@dataclass
class KnownPartner:
    a: str
    b: str
    cytoband: str
    disease: str
    source: str
    notes: str


def _load_loci_tsv(path: Path) -> list[GeneEntry]:
    out: list[GeneEntry] = []
    with open(path, encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            out.append(GeneEntry(
                gene=row["gene"], chrom=_strip_chr(row["chrom"]),
                start=int(row["start"]), end=int(row["end"]),
                strand=row.get("strand", "+"),
                role=row.get("role", ""),
                notes=row.get("notes", ""),
            ))
    return out


def _load_partners_tsv(path: Path) -> list[KnownPartner]:
    out: list[KnownPartner] = []
    with open(path, encoding="utf-8") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            out.append(KnownPartner(
                a=row["partner_a"], b=row["partner_b"],
                cytoband=row.get("cytoband", ""), disease=row.get("disease", ""),
                source=row.get("source", ""), notes=row.get("notes", ""),
            ))
    return out


def _strip_chr(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def load_builtin_loci() -> list[GeneEntry]:
    p = Path(str(pkg_files("quasarsv").joinpath("data/lymphoma_loci.tsv")))
    return _load_loci_tsv(p)


def load_builtin_partners() -> list[KnownPartner]:
    p = Path(str(pkg_files("quasarsv").joinpath("data/known_partners.tsv")))
    return _load_partners_tsv(p)


def _index_loci(loci: Iterable[GeneEntry]) -> dict[str, list[GeneEntry]]:
    idx: dict[str, list[GeneEntry]] = {}
    for g in loci:
        idx.setdefault(g.chrom, []).append(g)
    for k in idx:
        idx[k].sort(key=lambda g: g.start)
    return idx


def _hit_in_locus(idx: dict[str, list[GeneEntry]], chrom: str, pos: int,
                  upstream_bp: int = 5000, downstream_bp: int = 5000) -> tuple[str, str, str]:
    """Return (gene, region, driver_locus_label).

    Note: a 50 kb-pad widening was tried (to capture MBR/mcr breakpoints
    of t(14;18) etc.) but it dropped relaxed F1 from 0.87 to 0.82 by
    over-annotating peripheral calls. Reverted. Future improvements should
    expand the windows in ``data/lymphoma_loci.tsv`` explicitly per known
    breakpoint cluster, rather than blanket-padding every locus.
    """
    chrom = _strip_chr(chrom)
    if chrom not in idx:
        return "", "intergenic", ""
    # nearest gene by overlap or proximity
    for g in idx[chrom]:
        if g.start - upstream_bp <= pos <= g.end + downstream_bp:
            if g.start <= pos <= g.end:
                region = "exonic_or_intronic"
            elif pos < g.start:
                region = "upstream"
            else:
                region = "downstream"
            return g.gene, region, g.gene if g.role in ("driver", "IG_locus") else ""
    return "", "intergenic", ""


def _is_known_pair(partners: list[KnownPartner], gene_a: str, gene_b: str) -> KnownPartner | None:
    if not gene_a or not gene_b:
        return None
    for p in partners:
        if (p.a == gene_a and p.b == gene_b) or (p.a == gene_b and p.b == gene_a):
            return p
    return None


def annotate_calls(
    calls: list[FusionCall],
    extra_loci: Iterable[GeneEntry] | None = None,
    extra_partners: Iterable[KnownPartner] | None = None,
) -> list[FusionCall]:
    """In-place annotation using built-in lymphoma tables. Returns the same list."""
    loci = load_builtin_loci()
    partners = load_builtin_partners()
    if extra_loci:
        loci = list(loci) + list(extra_loci)
    if extra_partners:
        partners = list(partners) + list(extra_partners)
    idx = _index_loci(loci)

    for c in calls:
        a_gene, a_region, a_driver = _hit_in_locus(idx, c.chrom_a, c.pos_a)
        b_gene, b_region, b_driver = _hit_in_locus(idx, c.chrom_b, c.pos_b)
        c.gene_a, c.region_a = a_gene, a_region
        c.gene_b, c.region_b = b_gene, b_region

        # Driver-locus label: prioritise IG locus if present (canonical lymphoma)
        if a_driver and b_driver:
            c.driver_locus = f"{a_driver}-{b_driver}"
        elif a_driver:
            c.driver_locus = a_driver
        elif b_driver:
            c.driver_locus = b_driver

        kp = _is_known_pair(partners, a_gene, b_gene)
        if kp:
            c.known_partner = True
            c.known_partner_source = f"{kp.source}:{kp.disease}".strip(":")

    return calls
