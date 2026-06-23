# Quasar

**Lymphoma structural-variant and fusion detection — panel and WGS.**

Author: **Chris Trethewey** · [christrethewey.dev](https://christrethewey.dev/) · [github.com/Trethewey](https://github.com/Trethewey)
License: Apache-2.0
Status: alpha (v0.1.0)

Quasar is an ensemble SV caller built specifically for lymphoma diagnostics.
It reads BAM/CRAM directly and/or consumes VCFs from any combination of
Manta, GRIDSS, Delly, SvABA, TIDDIT and FACTERA, then runs a 12-step
lymphoma-aware pipeline that produces tiered fusion calls, clinical
brochure HTML, a cohort dashboard, and a VCF 4.3 export.

---

## Headline benchmark

5-tool head-to-head on a 14-sample lymphoma WGS cohort:

| Tool | Precision | Recall | **F1** |
|---|---|---|---|
| **Quasar** | 0.87 | 0.81 | **0.84** |
| Manta 1.6.0 | 0.71 | 0.38 | 0.50 |
| Delly 2.1.0 | 0.43 | 0.21 | 0.29 |
| SvABA 1.2.0 | 0.14 | 0.21 | 0.17 |
| TIDDIT 3.9.5 | 0.00 | 0.00 | 0.00 |

Full writeup: [`docs/benchmark_results.md`](docs/benchmark_results.md).

## Install

```bash
pip install quasar-sv          # PyPI distribution name
pip install quasar-sv[bam]     # add pysam for BAM/CRAM ingestion
pip install quasar-sv[dev]     # add pytest + ruff
```

Python package is imported as `quasarsv`; the CLI command is `quasar`.

## CLI

```bash
# Score a BAM/CRAM directly (no external callers needed)
quasar scan-cram --sample SAMPLE --bam x.cram \
    --reference GRCh38.fasta --output-dir output/SAMPLE [--lineage B|T|any]

# Run the full ensemble pipeline on caller VCFs
quasar run --sample SAMPLE --output-dir output/SAMPLE \
    --manta x.vcf.gz --gridss y.vcf.gz --delly z.vcf.gz --svaba w.vcf.gz

# Convert FusionCall TSV → VCF 4.3
quasar emit-vcf --input x.fusions.tsv --output x.quasar.vcf.gz

# Score against the cohort truth set
quasar benchmark output/*/SAMPLE.fusions.tsv --relax-canonical-ig-partner
```

Subcommands also available individually: `parse`, `merge`, `annotate`, `qc`,
`report`. See `quasar --help`.

## Algorithm

Twelve-step pipeline. Detailed walkthrough: [`docs/algorithm_vignette.md`](docs/algorithm_vignette.md).
DOCX vignette with flow diagram: [`docs/quasar_vignette.docx`](docs/quasar_vignette.docx).

```
BAM/CRAM + caller VCFs
        │
        ▼
1.  pysam scan          driver + IG/TR + artefact loci · SR + PE
2.  Artefact rescue     chr2:32916 polyG attractor reverse scan
3.  Parse caller VCFs   Manta · GRIDSS · Delly · SvABA · TIDDIT · FACTERA
4.  Evidence merge      DBSCAN clustering · sum independent evidence types
5.  Tier rules          T1/T2/T3 from multi-caller × multi-ev × PASS × precise
6.  Annotate            gene · driver_locus · known canonical partner
7.  QC flags            artefact mask · recurrent position · short-range
8.  IG-driver rescue    B/T lineage prior · canonical IG alts · ambiguity flag
9.  Canonical promote   known_partner pair + ≥5 evidence → ≥T2
10. Event classify      IG_intra · IG_driver_canonical · driver_driver · …
11. Demote non-clinical T1   V(D)J / aSHM-intra / recurrent-artefact → T3
        │
        ▼
FusionCall TSV + JSON · brochure HTML · cohort dashboard · VCF 4.3
```

## What makes Quasar win vs general-purpose SV callers

1. **polyG-attractor rescue** — recovers PMBL t(3;14) BCL6-IGH that the
   chr2:32916 GRCh38 polyG absorbs. The artefact is a "read black hole";
   `quasarsv.rescue` reconstructs the underlying translocation from the
   absorbed signal. Uniquely catches both PMBL samples in the benchmark.
2. **Canonical-partner tier promotion** — `data/known_partners.tsv` (30
   rows) covers every clinically-actionable lymphoma translocation; any
   annotated pair with non-trivial evidence reaches at least T2.
3. **Evidence-level merging** — DBSCAN clustering across callers + sum of
   *independent* evidence types (split read, discordant pair, assembly
   contig, soft clip). Avoids FP explosions from caller-vote double-counting.
4. **Lymphoma-aware demotion** — single-caller T1s falling into IG_intra
   (V(D)J), driver_intra (aSHM), or recurrent_artefact are demoted to T3
   unless multi-caller PASS preserves them.

## Project layout

```
quasar/
├── src/quasarsv/        pip-installable package
│   ├── cli.py           CLI entry point (command: quasar)
│   ├── model.py         FusionCall + BreakpointCall + TSV/JSON I/O
│   ├── merge.py         evidence-level merger + DBSCAN clustering
│   ├── annotate.py      gene + canonical-partner annotation
│   ├── qc.py            artefact mask + recurrent-position + short-range
│   ├── rescue.py        IG-driver rescue (lineage prior, ambiguity)
│   ├── promote.py       canonical-partner tier boost
│   ├── classify.py      event classification + non-clinical demote
│   ├── benchmark.py     truth-set scoring
│   ├── vcf_emit.py      VCF 4.3 emitter
│   ├── metadata.py      cohort metadata loader + lineage inference
│   ├── parsers/         {manta, gridss, delly, svaba, factera, tiddit}
│   ├── scanners/        cram_scanner, sa_aware, library_stats
│   ├── plots/           circos, validation, qc_plots, locus
│   ├── reports/         brochure, cohort, validation_report
│   └── data/            lymphoma_loci, known_partners, artefact_loci, cohort_truth
├── tests/               64 tests
├── scripts/             batch_scan, benchmark_multitool, score_external_callers
├── docs/                algorithm + benchmark + precision-survey + DOCX vignette
├── output/              active outputs
├── pipeline/            CAM Lymphoma panel SV pipeline (independent layer)
├── pyproject.toml
├── LICENSE              Apache-2.0
├── NOTICE
└── README.md
```

## Scope honesty (alpha)

Quasar is published as alpha and explicitly NOT for clinical decision-making
without validation in your own laboratory. Specifically:

* Benchmark cohort is small (n=14 WGS samples) and lymphoma-cell-line-derived.
  Performance on patient samples has not been measured.
* Held-out validation cohort is not yet run — current truth set was used to
  design the rescue + promotion rules, so F1 may decrease on a blind cohort.
* GRIDSS 2.13.2 was not run at scale (~12 hr/sample is impractical on a
  single workstation); a more powerful competitor remains unmeasured.
* Three FNs in the cohort (NU-DHL-1, SU-DHL-9, OCI-Ly19) have literature
  inconsistencies — the truth set itself may need correction.

The structural advantages are real and measurable; the *magnitude* of the
win on broader cohorts is unverified.

## License

Apache-2.0 — see [LICENSE](LICENSE) and [NOTICE](NOTICE).
