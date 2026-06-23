#!/usr/bin/env bash
# wgs_sv_callers.sh — produce SV VCFs from a WGS CRAM for quasarsv ingest.
#
# Run from WSL or Linux. Tools required on PATH:
#   - samtools (>= 1.17)
#   - manta (>= 1.6)
#   - configManta.py / runWorkflow.py in PATH
#   - gridss (>= 2.13)
#   - delly (>= 1.2)
#   - svaba (>= 1.1)
#
# Usage:
#   wgs_sv_callers.sh \
#       --sample <sample_id> \
#       --cram   <path/to/sample.cram> \
#       --ref    <path/to/GRCh38_full_with_decoys.fa> \
#       --out    <path/to/output/dir> \
#       [--threads 16] [--callers manta,gridss,delly,svaba]
#
# Produces, for each caller, $OUT/$SAMPLE.<caller>.vcf.gz suitable for:
#   quasarsv run --sample $SAMPLE --output-dir $OUT \
#                   --manta  $OUT/$SAMPLE.manta.vcf.gz \
#                   --gridss $OUT/$SAMPLE.gridss.vcf.gz \
#                   --delly  $OUT/$SAMPLE.delly.vcf.gz \
#                   --svaba  $OUT/$SAMPLE.svaba.vcf.gz
#
# Note on reference: the WGS CRAMs in W:\WGS_data use chr-prefixed contigs and
# the full GRCh38 decoy/ALT set (3,366 contigs). They are NOT compatible with
# the panel reference. Use the matched reference for the source cohort.

set -euo pipefail

THREADS=16
CALLERS="manta,gridss,delly,svaba"
SAMPLE=""; CRAM=""; REF=""; OUT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample)  SAMPLE="$2"; shift 2;;
    --cram)    CRAM="$2";   shift 2;;
    --ref)     REF="$2";    shift 2;;
    --out)     OUT="$2";    shift 2;;
    --threads) THREADS="$2";shift 2;;
    --callers) CALLERS="$2";shift 2;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *)         echo "unknown arg $1" >&2; exit 2;;
  esac
done

for v in SAMPLE CRAM REF OUT; do
  [[ -n "${!v}" ]] || { echo "missing --${v,,}" >&2; exit 2; }
done

mkdir -p "$OUT" "$OUT/logs"
log() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

run_manta() {
  log "manta: configure"
  rm -rf "$OUT/manta_run"
  configManta.py --bam "$CRAM" --referenceFasta "$REF" --runDir "$OUT/manta_run" \
    >> "$OUT/logs/manta.log" 2>&1
  log "manta: workflow"
  "$OUT/manta_run/runWorkflow.py" -j "$THREADS" \
    >> "$OUT/logs/manta.log" 2>&1
  cp "$OUT/manta_run/results/variants/candidateSV.vcf.gz" \
     "$OUT/${SAMPLE}.manta.vcf.gz"
  cp "$OUT/manta_run/results/variants/candidateSV.vcf.gz.tbi" \
     "$OUT/${SAMPLE}.manta.vcf.gz.tbi"
  log "manta: done"
}

run_gridss() {
  log "gridss: starting"
  gridss -t "$THREADS" -r "$REF" \
    -o "$OUT/${SAMPLE}.gridss.vcf.gz" \
    -a "$OUT/${SAMPLE}.gridss.assembly.bam" \
    --workingdir "$OUT/gridss_work" \
    "$CRAM" \
    >> "$OUT/logs/gridss.log" 2>&1
  log "gridss: done"
}

run_delly() {
  log "delly: starting"
  delly call -g "$REF" -o "$OUT/${SAMPLE}.delly.bcf" "$CRAM" \
    >> "$OUT/logs/delly.log" 2>&1
  bcftools view "$OUT/${SAMPLE}.delly.bcf" | bgzip > "$OUT/${SAMPLE}.delly.vcf.gz"
  tabix -p vcf "$OUT/${SAMPLE}.delly.vcf.gz"
  log "delly: done"
}

run_svaba() {
  log "svaba: starting"
  svaba run -t "$CRAM" -G "$REF" -p "$THREADS" \
    -a "$OUT/${SAMPLE}.svaba" -L 6 \
    >> "$OUT/logs/svaba.log" 2>&1
  # SvABA emits *.svaba.sv.vcf and small variants .indel.vcf — fusion path uses the SV one
  bgzip -c "$OUT/${SAMPLE}.svaba.svaba.sv.vcf" > "$OUT/${SAMPLE}.svaba.vcf.gz"
  tabix -p vcf "$OUT/${SAMPLE}.svaba.vcf.gz"
  log "svaba: done"
}

IFS=',' read -ra CARR <<< "$CALLERS"
for c in "${CARR[@]}"; do
  case "$c" in
    manta)  run_manta  ;;
    gridss) run_gridss ;;
    delly)  run_delly  ;;
    svaba)  run_svaba  ;;
    *) echo "unknown caller $c" >&2; exit 3;;
  esac
done

log "all done. Next: run quasarsv with the produced VCFs."
echo
echo "    quasarsv run --sample $SAMPLE --output-dir $OUT \\"
echo "        --manta  $OUT/${SAMPLE}.manta.vcf.gz \\"
echo "        --gridss $OUT/${SAMPLE}.gridss.vcf.gz \\"
echo "        --delly  $OUT/${SAMPLE}.delly.vcf.gz \\"
echo "        --svaba  $OUT/${SAMPLE}.svaba.vcf.gz"
