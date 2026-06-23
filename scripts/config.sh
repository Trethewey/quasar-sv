#!/usr/bin/env bash
# Central configuration for quasarsv scripts.
# Source this from the top of any script:  source "$(dirname "$0")/config.sh"
#
# Every path here can be overridden via environment variable. Defaults are
# discoverable via `$PROJECT_ROOT`-relative lookups.

# Project root is two directories up from this file (scripts/config.sh).
# Works whether the project lives at /mnt/d/..., /home/user/..., or elsewhere.
PROJECT_ROOT="${PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# Reference genome (must be indexed: .fai, .dict, .bwt for SvABA/GRIDSS).
# Override with FUSIONFORGE_REFERENCE=/your/path/Homo_sapiens_assembly38.fasta
FF_REFERENCE="${FUSIONFORGE_REFERENCE:-$PROJECT_ROOT/ref/Homo_sapiens_assembly38.fasta}"

# CRAM source root — where the cohort CRAMs live (typically a mounted NAS).
# Override with FUSIONFORGE_CRAM_ROOT=/your/path/WGS_data
FF_CRAM_ROOT="${FUSIONFORGE_CRAM_ROOT:-$PROJECT_ROOT/cram_source}"

# Local CRAM cache — fast-local copy used by the benchmark harness to
# bypass slow filesystem layers. Override with FUSIONFORGE_LOCAL_CRAM_DIR.
FF_LOCAL_CRAM_DIR="${FUSIONFORGE_LOCAL_CRAM_DIR:-$PROJECT_ROOT/cram_local}"

# Cohort metadata XLSX with one row per sample. Optional; lineage inference
# falls back to "B" if absent. Override with FUSIONFORGE_COHORT_METADATA.
FF_COHORT_METADATA="${FUSIONFORGE_COHORT_METADATA:-$FF_CRAM_ROOT/cohort_metadata.xlsx}"

# Conda / mamba installation root. Override with FUSIONFORGE_CONDA_ROOT.
FF_CONDA_ROOT="${FUSIONFORGE_CONDA_ROOT:-$HOME/miniconda3}"

# Thread count for tools that accept it.
FF_THREADS="${FUSIONFORGE_THREADS:-8}"

# Where benchmark + cohort outputs land.
FF_OUTPUT_DIR="${FUSIONFORGE_OUTPUT_DIR:-$PROJECT_ROOT/output}"

# Lymphoma target regions BED (driver + IG/TR + artefact loci, ±200 kb pad).
FF_TARGETS_BED="${FUSIONFORGE_TARGETS_BED:-$FF_OUTPUT_DIR/config/lymphoma_targets.bed}"

# Python path so `python3 -m quasarsv.cli ...` works without install.
export PYTHONPATH="${PROJECT_ROOT}/src:${PYTHONPATH:-}"

# Export everything so child processes see it.
export PROJECT_ROOT FF_REFERENCE FF_CRAM_ROOT FF_LOCAL_CRAM_DIR \
       FF_COHORT_METADATA FF_CONDA_ROOT FF_THREADS FF_OUTPUT_DIR \
       FF_TARGETS_BED
