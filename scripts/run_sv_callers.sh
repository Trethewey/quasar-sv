#!/usr/bin/env bash
# run_sv_callers.sh — run each external SV caller on one CRAM, write per-tool VCFs.
#
# Tools driven (one mamba env each, all from bioconda):
#   manta   (Manta 1.6.0)
#   delly   (Delly 2.1.0)
#   svaba   (SvABA 1.2.0)
#   gridss  (GRIDSS 2.13)
#   tiddit  (TIDDIT 3.9.5)
#   lumpy   (Lumpy-SV)
#
# Outputs:
#   $OUT/manta/$SAMPLE/results/variants/diploidSV.vcf.gz
#   $OUT/delly/$SAMPLE/$SAMPLE.delly.bcf
#   $OUT/svaba/$SAMPLE/$SAMPLE.svaba.svaba.sv.vcf
#   $OUT/gridss/$SAMPLE/$SAMPLE.gridss.vcf.gz
#   $OUT/tiddit/$SAMPLE/$SAMPLE.tiddit.vcf
#   $OUT/lumpy/$SAMPLE/$SAMPLE.lumpy.vcf
#
# Resume-aware: skips any caller whose final VCF already exists.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

PROJECT="$PROJECT_ROOT"
THREADS=${THREADS:-$FF_THREADS}

SAMPLE=""
CRAM=""
REF=""
OUT=""
TOOLS="manta,delly,svaba,gridss,tiddit"
LOG_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)  SAMPLE="$2"; shift 2;;
    --cram)    CRAM="$2"; shift 2;;
    --ref)     REF="$2"; shift 2;;
    --out)     OUT="$2"; shift 2;;
    --tools)   TOOLS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) sed -n '1,30p' "$0"; exit 0;;
    *) echo "unknown arg $1" >&2; exit 2;;
  esac
done

for v in SAMPLE CRAM REF OUT; do
  [[ -n "${!v}" ]] || { echo "missing --${v,,}" >&2; exit 2; }
done

[[ -f "$CRAM" ]] || { echo "CRAM not found: $CRAM" >&2; exit 3; }
[[ -f "$REF" ]] || { echo "REF not found: $REF" >&2; exit 3; }

LOG_DIR=$OUT/logs
mkdir -p "$LOG_DIR"
log() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

export PATH="$HOME/miniconda3/bin:$PATH"

run_manta() {
  local dir=$OUT/manta/$SAMPLE
  local vcf=$dir/results/variants/diploidSV.vcf.gz
  if [[ -s "$vcf" ]]; then
    log "manta: $SAMPLE — already done"
    return 0
  fi
  log "manta: $SAMPLE — configure + run"
  mkdir -p "$dir"
  rm -rf "$dir/runDir"
  mamba run -n manta configManta.py \
    --bam "$CRAM" \
    --referenceFasta "$REF" \
    --runDir "$dir/runDir" >"$LOG_DIR/manta_${SAMPLE}_configure.log" 2>&1
  mamba run -n manta python2 "$dir/runDir/runWorkflow.py" \
    -j "$THREADS" -m local >"$LOG_DIR/manta_${SAMPLE}_run.log" 2>&1
  mkdir -p "$dir/results/variants"
  cp "$dir/runDir/results/variants/diploidSV.vcf.gz"     "$dir/results/variants/" 2>/dev/null || true
  cp "$dir/runDir/results/variants/diploidSV.vcf.gz.tbi" "$dir/results/variants/" 2>/dev/null || true
  log "manta: $SAMPLE — done"
}

run_delly() {
  local dir=$OUT/delly/$SAMPLE
  local bcf=$dir/$SAMPLE.delly.bcf
  if [[ -s "$bcf" ]]; then
    log "delly: $SAMPLE — already done"
    return 0
  fi
  log "delly: $SAMPLE — calling all SV types"
  mkdir -p "$dir"
  OMP_NUM_THREADS=$THREADS mamba run -n delly delly call \
    -g "$REF" \
    -o "$bcf" \
    "$CRAM" >"$LOG_DIR/delly_${SAMPLE}.log" 2>&1
  # Convert BCF -> VCF for quasarsv ingest
  mamba run -n base-bio bcftools view "$bcf" \
    -O z -o "$dir/$SAMPLE.delly.vcf.gz" 2>"$LOG_DIR/delly_${SAMPLE}_bcf.log" || true
  mamba run -n base-bio bcftools index -t "$dir/$SAMPLE.delly.vcf.gz" 2>/dev/null || true
  log "delly: $SAMPLE — done"
}

run_svaba() {
  local dir=$OUT/svaba/$SAMPLE
  local vcf=$dir/$SAMPLE.svaba.svaba.sv.vcf
  if [[ -s "$vcf" ]]; then
    log "svaba: $SAMPLE — already done"
    return 0
  fi
  log "svaba: $SAMPLE — local assembly + SV calling"
  mkdir -p "$dir"
  pushd "$dir" >/dev/null
  mamba run -n svaba svaba run \
    -t "$CRAM" \
    -G "$REF" \
    -a "$SAMPLE.svaba" \
    -p "$THREADS" >"$LOG_DIR/svaba_${SAMPLE}.log" 2>&1
  popd >/dev/null
  log "svaba: $SAMPLE — done"
}

run_gridss() {
  local dir=$OUT/gridss/$SAMPLE
  local vcf=$dir/$SAMPLE.gridss.vcf.gz
  if [[ -s "$vcf" ]]; then
    log "gridss: $SAMPLE — already done"
    return 0
  fi
  log "gridss: $SAMPLE — assembly + SV calling (slowest)"
  mkdir -p "$dir"
  # GRIDSS needs an indexed BAM (not CRAM) and bwa-indexed reference
  # Convert CRAM -> BAM (random-access; only emits what's needed)
  local bam=$dir/$SAMPLE.bam
  if [[ ! -s "$bam" ]]; then
    log "gridss: converting CRAM -> BAM (needed once)"
    mamba run -n base-bio samtools view -b -@ "$THREADS" -T "$REF" -o "$bam" "$CRAM"
    mamba run -n base-bio samtools index -@ "$THREADS" "$bam"
  fi
  mamba run -n gridss gridss \
    -r "$REF" \
    -o "$vcf" \
    -a "$dir/$SAMPLE.assembly.bam" \
    -t "$THREADS" \
    --workingdir "$dir/workingdir" \
    "$bam" >"$LOG_DIR/gridss_${SAMPLE}.log" 2>&1
  log "gridss: $SAMPLE — done"
}

run_tiddit() {
  local dir=$OUT/tiddit/$SAMPLE
  local vcf=$dir/$SAMPLE.tiddit.vcf
  if [[ -s "$vcf" ]]; then
    log "tiddit: $SAMPLE — already done"
    return 0
  fi
  log "tiddit: $SAMPLE — coverage-based SV calling"
  mkdir -p "$dir"
  pushd "$dir" >/dev/null
  mamba run -n tiddit tiddit --sv \
    --bam "$CRAM" \
    --ref "$REF" \
    -o "$SAMPLE.tiddit" \
    --threads "$THREADS" >"$LOG_DIR/tiddit_${SAMPLE}.log" 2>&1
  [[ -f "$SAMPLE.tiddit.vcf" ]] && mv "$SAMPLE.tiddit.vcf" "$vcf" || true
  popd >/dev/null
  log "tiddit: $SAMPLE — done"
}

run_lumpy() {
  local dir=$OUT/lumpy/$SAMPLE
  local vcf=$dir/$SAMPLE.lumpy.vcf
  if [[ -s "$vcf" ]]; then
    log "lumpy: $SAMPLE — already done"
    return 0
  fi
  log "lumpy: $SAMPLE — extracting discordant + split + calling"
  mkdir -p "$dir"
  # Lumpy needs split/discordant BAMs prepared
  local bam=$OUT/gridss/$SAMPLE/$SAMPLE.bam
  if [[ ! -s "$bam" ]]; then
    log "lumpy: BAM dependency missing — run gridss first to materialise BAM"
    return 1
  fi
  local discordant=$dir/$SAMPLE.discordant.bam
  local split=$dir/$SAMPLE.split.bam
  mamba run -n base-bio samtools view -b -F 1294 "$bam" \
    | mamba run -n base-bio samtools sort -@ "$THREADS" -o "$discordant" -
  mamba run -n base-bio samtools index "$discordant"
  mamba run -n lumpy lumpyexpress \
    -B "$bam" \
    -S "$split" \
    -D "$discordant" \
    -o "$vcf" >"$LOG_DIR/lumpy_${SAMPLE}.log" 2>&1 || true
  log "lumpy: $SAMPLE — done (best-effort)"
}

IFS=',' read -ra TOOL_LIST <<< "$TOOLS"
for tool in "${TOOL_LIST[@]}"; do
  case "$tool" in
    manta)  run_manta;;
    delly)  run_delly;;
    svaba)  run_svaba;;
    gridss) run_gridss;;
    tiddit) run_tiddit;;
    lumpy)  run_lumpy;;
    *) echo "unknown tool: $tool" >&2;;
  esac
done

log "ALL DONE: $SAMPLE"
