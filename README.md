<p align="center">
  <img src="brand/png/quasar-banner-dark.png#gh-dark-mode-only" alt="Quasar" width="720">
  <img src="brand/png/quasar-banner-light.png#gh-light-mode-only" alt="Quasar" width="720">
</p>

<br>

<p align="center">
  <a href="LICENSE"><img alt="Licence" src="https://img.shields.io/badge/licence-Apache_2.0-blue.svg"></a>
  <a href="pyproject.toml"><img alt="Python" src="https://img.shields.io/badge/python-3.10%2B-blue.svg"></a>
  <a href="#scope-and-limitations"><img alt="Status" src="https://img.shields.io/badge/status-alpha-orange.svg"></a>
</p>

Quasar is a structural-variant caller for next-generation sequencing data,
specialised for lymphoma diagnostics. It identifies clinically actionable
translocations — including IGH-BCL2, IGH-MYC, BCL6-IGH, NPM1-ALK and the
canonical lymphoma partner panel — from aligned reads in BAM or CRAM,
and emits results as VCF 4.3 alongside a per-sample clinical report.

## Overview

Quasar combines a pysam-based read scanner with lymphoma-specific
annotation and tier-promotion logic. Three design choices distinguish
it from general-purpose SV callers:

**Artefact-aware translocation rescue.** Chimeric reads from t(3;14)
BCL6-IGH and related IG-switch events route into a known GRCh38
reference artefact at chr2:32916, where they are absorbed by a
poly-G motif and lost to ordinary detection. Quasar's rescue layer
reconstructs the underlying translocation from artefact-mediated
signal.

**Canonical-partner annotation.** A curated table of lymphoma
translocation partners drives clinical-tier promotion of annotated
pairs, including events that general-purpose callers filter as
low-quality noise.

**Non-clinical event demotion.** Recurrent V(D)J recombination events
and aSHM-hotspot signals are demoted from the clinical tier unless
corroborated by multiple independent callers.

## Installation

```bash
pip install quasar-sv
```

CLI entry point: `quasar`. Python package: `quasarsv`.

## Quick start

```bash
quasar call \
    --sample SAMPLE_ID \
    --bam aligned.cram \
    --reference GRCh38.fasta \
    --output-dir output/SAMPLE_ID
```

## Output

Files written to `output/SAMPLE_ID/`:

| File | Description |
|------|-------------|
| `SAMPLE_ID.fusions.tsv` | All calls with tier, partner genes, evidence summary, QC flags |
| `SAMPLE_ID.fusions.vcf.gz` | Calls in VCF 4.3 format |
| `brochure_SAMPLE_ID.html` | Per-sample clinical brochure |

Calls are stratified into three confidence tiers:

| Tier | Definition |
|------|------------|
| **T1** | High confidence. Canonical lymphoma partner pair with strong supporting evidence, or independent PASS from multiple callers. |
| **T2** | Moderate confidence. Annotated lymphoma partner with reduced evidence, or PASS in two or more callers. |
| **T3** | Low confidence. Surfaced for review. |

## Documentation

| Document | Content |
|----------|---------|
| [`docs/algorithm_vignette.md`](docs/algorithm_vignette.md) | Full algorithm walkthrough |
| [`docs/benchmark_results.md`](docs/benchmark_results.md) | Head-to-head benchmark and caveats |
| [`docs/panel_validation.md`](docs/panel_validation.md) | Targeted panel BAM usage |
| [`docs/quasar_vignette.docx`](docs/quasar_vignette.docx) | Printable algorithm summary |

## Scope and limitations

Quasar v0.1.0 is an alpha release. Validation has been performed on a
14-sample lymphoma cell-line WGS cohort. Held-out cohorts, broad
cohort generalisability, and patient-sample performance have not been
assessed. Quasar is not intended for clinical decision-making without
site-specific validation.

## Citation

A manuscript is in preparation. In the interim, please cite this
repository:

```
Trethewey C. Quasar: structural-variant calling for lymphoma sequencing.
https://github.com/Trethewey/Quasar (2026).
```

## Author

Chris Trethewey · [christrethewey.dev](https://christrethewey.dev/) · [github.com/Trethewey](https://github.com/Trethewey)

## Licence

Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
