<p align="center">
  <img src="brand/quasar-logo-dark-stacked.svg" alt="Quasar" width="520">
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

**Artefact-aware breakend handling.** A known GRCh38 poly-G attractor at
chr2:32,916 absorbs adapter read-through and 2-colour poly-G tails from
every locus at a near-constant rate, carrying no information about what
joins to what. Quasar filters these clips at source and marks any call
whose partner breakend falls in a masked artefact region, so an
unresolvable breakend is reported as `partner_undetermined` rather than
being assigned a partner it cannot measure.

**Canonical-partner annotation.** A curated table of lymphoma
translocation partners drives clinical-tier promotion of annotated
pairs, including events that general-purpose callers filter as
low-quality noise.

**Non-clinical event demotion.** Recurrent V(D)J recombination events
and aSHM-hotspot signals are demoted from the clinical tier unless
corroborated by multiple independent callers.

## Installation

From a clone (the repository is currently private, so this is the reliable route):

```bash
git clone https://github.com/Trethewey/Quasar.git
cd Quasar
pip install ".[bam]"
```

Or directly, if you have repository access:

```bash
pip install "quasar-sv[bam] @ git+ssh://git@github.com/Trethewey/Quasar.git"
```

The `[bam]` extra pulls in pysam, without which the BAM/CRAM scanner cannot run.

CLI entry point: `quasar`. Python package: `quasarsv`. (A PyPI release under the
name `quasar-sv` is planned; until then install from the repository as above.)

**Requirements.** Quasar reads BAM/CRAM via [pysam](https://github.com/pysam-developers/pysam),
which ships wheels for **Linux and macOS only** — on Windows, run under WSL.
The reference FASTA must be indexed (`samtools faidx GRCh38.fasta`, producing
`GRCh38.fasta.fai`) so the scanner can fetch loci by region.

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

No head-to-head benchmark is published. An earlier one was withdrawn: every
external caller in it had been run on region-restricted input by our own
harness, so their scores measured our configuration rather than their
performance. Fair genome-wide re-runs are in progress; no comparative claim
will be made until they are complete and scored under identical rules.

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
