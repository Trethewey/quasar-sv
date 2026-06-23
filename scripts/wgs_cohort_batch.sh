#!/usr/bin/env bash
# wgs_cohort_batch.sh — drive wgs_sv_callers.sh across the W:\WGS_data cohort,
# pulling sample identities from cohort_metadata.xlsx via a flat samples.tsv.
#
# Expected samples.tsv columns (tab-separated, no header):
#   sample_id   cram_path   reference_fa
#
# Each line is dispatched serially (callers within a line use --threads).
# Logs go to $OUT_ROOT/<sample>/logs/.

set -euo pipefail

SAMPLES_TSV=""; OUT_ROOT=""; THREADS=16; CALLERS="manta,gridss,delly,svaba"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples)  SAMPLES_TSV="$2"; shift 2;;
    --out-root) OUT_ROOT="$2";    shift 2;;
    --threads)  THREADS="$2";     shift 2;;
    --callers)  CALLERS="$2";     shift 2;;
    -h|--help)  sed -n '2,12p' "$0"; exit 0;;
    *)          echo "unknown arg $1" >&2; exit 2;;
  esac
done

[[ -n "$SAMPLES_TSV" && -n "$OUT_ROOT" ]] || {
  echo "usage: $0 --samples samples.tsv --out-root /path/out [--threads N]" >&2
  exit 2
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while IFS=$'\t' read -r sample cram ref; do
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue
  out="$OUT_ROOT/$sample"
  echo "==== $sample → $out"
  "$SCRIPT_DIR/wgs_sv_callers.sh" \
    --sample "$sample" --cram "$cram" --ref "$ref" --out "$out" \
    --threads "$THREADS" --callers "$CALLERS"

  # Run quasarsv immediately after callers finish so each sample's report
  # lands in its own dir (parallelisable across samples — wrap with GNU parallel).
  quasarsv run --sample "$sample" --output-dir "$out" \
    --manta  "$out/${sample}.manta.vcf.gz" \
    --gridss "$out/${sample}.gridss.vcf.gz" \
    --delly  "$out/${sample}.delly.vcf.gz" \
    --svaba  "$out/${sample}.svaba.vcf.gz"
done < "$SAMPLES_TSV"

echo "cohort complete. Render cross-sample dashboards with:"
echo "    quasarsv report --input $OUT_ROOT/cohort.fusions.tsv --output-dir $OUT_ROOT/cohort --kind cohort"
