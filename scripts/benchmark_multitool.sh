#!/usr/bin/env bash
# Run external SV callers on the cohort and score each against the truth set.
# Resume-aware. Override via env: SAMPLES, TOOLS, TARGETS, THREADS, CRAM_DIR.

set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

PROJECT="$PROJECT_ROOT"
REF="$FF_REFERENCE"
OUT="$FF_OUTPUT_DIR/benchmark"
COHORT_OUT="$FF_OUTPUT_DIR/wgs_cohort"
TARGETS=${TARGETS:-$FF_TARGETS_BED}
THREADS=${THREADS:-$FF_THREADS}
TOOLS=${TOOLS:-manta,delly,svaba,tiddit}
TRUTH=${TRUTH:-$PROJECT/src/quasarsv/data/cohort_truth.tsv}

mkdir -p "$OUT/logs"
LOG=$OUT/_run.log
log() {
  # IMPORTANT: write to stderr (not stdout) so that helper functions which
  # capture their result via `$(...)` aren't polluted by progress messages.
  # E.g. _prepare_targeted_bam's `bam=$(_prepare_targeted_bam ...)` would
  # otherwise include every log line in the returned string.
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >> "$LOG"
  printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*" >&2
}

# Default positive-truth samples (skip negative controls + unscanned)
DEFAULT_SAMPLES="ERR9128954_U2940 ERR9188549_Karpas1106P SRR1236466_OCI-Ly1_DLBCL_cell_line SRR1236467_NU-DHL-1_DLBCL_cell_line SRR1236468_DB_DLBCL_cell_line SRR1236469_NU-DUL-1_DLBCL_cell_line SRR1236470_OCI-Ly7_DLBCL_cell_line SRR1236472_MD903_DLBCL_cell_line SRR1236473_SU-DHL-9_DLBCL_cell_line SRR1236474_DOHH-2_DLBCL_cell_line SRR1236475_WSU-DLCL2_DLBCL_cell_line SRR1236476_SU-DHL-6_DLBCL_cell_line SRR1236477_OCI-Ly19_DLBCL_cell_line SRR1236478_Karpas-422_DLBCL_cell_line SRR16382373_Human_sample_from_Homo_sapiens SRR16382375_Human_sample_from_Homo_sapiens"
SAMPLES=${SAMPLES:-$DEFAULT_SAMPLES}

# Read manifest: sample <tab> cram
declare -A CRAM_OF
while IFS=$'\t' read -r sample cram _; do
  [[ -z "$sample" || "$sample" =~ ^# ]] && continue
  CRAM_OF[$sample]="$cram"
done < "$COHORT_OUT/_samples.tsv"

# Optional CRAM_DIR override — when CRAMs have been pre-localized to a fast
# filesystem (e.g. F:\Data\cram_local mounted at /mnt/f/Data/cram_local), point
# every sample at the local copy and skip the pre-extract step entirely.
if [[ -n "${CRAM_DIR:-}" ]]; then
  log "CRAM_DIR override active: looking for *.cram in $CRAM_DIR"
  for sample in "${!CRAM_OF[@]}"; do
    # Extract the SRR/ERR/DRR accession from the sample id (first underscored token)
    accession="${sample%%_*}"
    candidate="$CRAM_DIR/${accession}.cram"
    if [[ -s "$candidate" ]]; then
      CRAM_OF[$sample]="$candidate"
    fi
  done
fi

export PATH="$FF_CONDA_ROOT/bin:$PATH"

# ----- per-caller runners ------------------------------------------------

run_manta_targeted() {
  local sample=$1 cram=$2 dir=$OUT/manta/$sample
  local vcf=$dir/results/variants/diploidSV.vcf.gz
  if [[ -s "$vcf" ]]; then echo "$vcf"; return 0; fi
  mkdir -p "$dir"
  rm -rf "$dir/runDir"
  local bedgz=$dir/targets.bed.gz
  cp "$TARGETS" "$dir/targets.bed"
  mamba run -n base-bio bgzip -f "$dir/targets.bed"
  mamba run -n base-bio tabix -p bed "$bedgz"
  mamba run -n manta configManta.py \
    --bam "$cram" --referenceFasta "$REF" \
    --callRegions "$bedgz" --runDir "$dir/runDir" \
    >"$OUT/logs/manta_${sample}_cfg.log" 2>&1
  mamba run -n manta python2 "$dir/runDir/runWorkflow.py" \
    -j "$THREADS" -m local >"$OUT/logs/manta_${sample}_run.log" 2>&1
  mkdir -p "$dir/results/variants"
  cp "$dir/runDir/results/variants/diploidSV.vcf.gz"     "$dir/results/variants/"
  cp "$dir/runDir/results/variants/diploidSV.vcf.gz.tbi" "$dir/results/variants/" 2>/dev/null || true
  echo "$vcf"
}

_prepare_targeted_bam() {
  # Extract ONLY lymphoma target regions from CRAM into a small local BAM.
  # One NAS read per sample; subsequent tool calls run against the local BAM.
  # BAM is shared across delly + tiddit per-sample (and torn down at end).
  #
  # Strategy: localize-first. samtools' CRAM index-based seeks are slow
  # on remote / SMB / overlay filesystems. Sequentially copy the CRAM to
  # a fast local filesystem, then run samtools against the local copy.
  #   1. cp CRAM to /tmp/ff_cramcache/   (one big sequential read)
  #   2. samtools view -L bed on the local CRAM
  #   3. delete the local CRAM after the targeted BAM is built
  # Retries up to 3 times on transient I/O failures.
  local sample=$1 cram=$2
  local bam_dir=/tmp/ff_targeted
  local cram_dir=/tmp/ff_cramcache
  mkdir -p "$bam_dir" "$cram_dir"
  local bam=$bam_dir/${sample}.targeted.bam
  if [[ -s "$bam" && -s "${bam}.bai" ]]; then
    echo "$bam"; return 0
  fi
  local local_cram=$cram_dir/${sample}.cram
  local local_crai=${local_cram}.crai

  # Step 1: localize the CRAM (and its CRAI if available). Fast path:
  # if the CRAM already lives on a fast local filesystem (via FF_LOCAL_CRAM_DIR
  # or CRAM_DIR override), skip the localize copy and use it directly.
  local cram_path="$cram"
  case "$cram" in
    /mnt/f/*|/mnt/d/*)
      log "  CRAM already on local fs: $cram (skipping localize)"
      # Look for the CRAI next to the CRAM
      local cram_local_crai="${cram}.crai"
      if [[ -f "$cram_local_crai" ]]; then
        cram_path="$cram"
      else
        log "  warning: no .crai next to local CRAM; samtools will need one"
      fi
      ;;
    *)
      for attempt in 1 2 3; do
        if [[ ! -s "$local_cram" ]]; then
          log "  localize CRAM ($sample, attempt $attempt): $cram -> $local_cram"
          if ! cp "$cram" "$local_cram" 2>"$OUT/logs/localize_${sample}.log"; then
            log "  cp failed; backoff $((10 * attempt))s and retry"
            rm -f "$local_cram"
            sleep $((10 * attempt))
            continue
          fi
        fi
        # Try to grab the matching CRAI; if missing we'll let samtools build it.
        local nas_crai="${cram}.crai"
        if [[ -f "$nas_crai" && ! -s "$local_crai" ]]; then
          cp "$nas_crai" "$local_crai" 2>>"$OUT/logs/localize_${sample}.log" || true
        fi
        break
      done
      cram_path="$local_cram"
      ;;
  esac
  if [[ -z "$cram_path" || ! -s "$cram_path" ]]; then
    log "  localize: ALL 3 retries exhausted for $sample"
    return 1
  fi

  # Step 2: targeted extract from LOCAL CRAM. Fast random-access via CRAI.
  log "  targeted extract from $cram_path -> $bam"
  if ! mamba run -n base-bio samtools view -b \
      -@ "$THREADS" -T "$REF" -L "$TARGETS" \
      -o "$bam" "$cram_path" 2>"$OUT/logs/preextract_${sample}.log"; then
    log "  samtools view -L failed for $sample"
    rm -f "$bam"
    return 1
  fi
  mamba run -n base-bio samtools index -@ "$THREADS" "$bam" 2>/dev/null || true
  local sz=$(stat -c '%s' "$bam" 2>/dev/null || echo 0)
  log "  $sample targeted BAM ready (${sz} bytes)"
  echo "$bam"
  return 0
}

_cleanup_targeted_bam() {
  local sample=$1
  rm -f "/tmp/ff_targeted/${sample}.targeted.bam" "/tmp/ff_targeted/${sample}.targeted.bam.bai" 2>/dev/null
  # Also reclaim the localized CRAM — at 60-85 GB it's the dominant disk user.
  rm -f "/tmp/ff_cramcache/${sample}.cram" "/tmp/ff_cramcache/${sample}.cram.crai" 2>/dev/null
}

run_delly_targeted() {
  local sample=$1 cram=$2 dir=$OUT/delly/$sample
  local bcf=$dir/$sample.delly.bcf
  local vcf=$dir/$sample.delly.vcf.gz
  if [[ -s "$vcf" ]]; then echo "$vcf"; return 0; fi
  mkdir -p "$dir"
  # Use the pre-extracted targeted BAM (one NAS read shared across tools).
  local bam
  bam=$(_prepare_targeted_bam "$sample" "$cram")
  [[ -z "$bam" || ! -s "$bam" ]] && { log "  delly: no targeted BAM, skipping"; echo ""; return 1; }
  OMP_NUM_THREADS=$THREADS mamba run -n delly delly call \
    -g "$REF" \
    -o "$bcf" "$bam" >"$OUT/logs/delly_${sample}.log" 2>&1 || true
  if [[ -s "$bcf" ]]; then
    mamba run -n base-bio bcftools view "$bcf" -O z -o "$vcf" \
      2>"$OUT/logs/delly_${sample}_bcf.log" || true
    mamba run -n base-bio bcftools index -t "$vcf" 2>/dev/null || true
  fi
  echo "$vcf"
}

run_gridss_targeted() {
  local sample=$1 cram=$2 dir=$OUT/gridss/$sample
  local vcf=$dir/$sample.gridss.vcf.gz
  if [[ -s "$vcf" ]]; then echo "$vcf"; return 0; fi
  mkdir -p "$dir"
  # GRIDSS requires an indexed BAM. We convert on the WSL rootfs (more room
  # than D:) and DELETE after the GRIDSS run to stay below disk budget.
  local bam=/tmp/gridss_${sample}.bam
  log "gridss: $sample — CRAM->BAM convert to /tmp (cleaned after)"
  mamba run -n base-bio samtools view -b -@ "$THREADS" -T "$REF" -o "$bam" "$cram" 2>"$OUT/logs/gridss_${sample}_bamconv.log" || { rm -f "$bam"; return 1; }
  mamba run -n base-bio samtools index -@ "$THREADS" "$bam" 2>/dev/null || true
  log "gridss: $sample — running"
  mamba run -n gridss gridss \
    -r "$REF" \
    -o "$vcf" \
    -a "$dir/$sample.assembly.bam" \
    -t "$THREADS" \
    --workingdir "$dir/workingdir" \
    -l "$TARGETS" \
    "$bam" >"$OUT/logs/gridss_${sample}.log" 2>&1
  local rc=$?
  log "gridss: $sample — cleaning $bam (regardless of run rc=$rc)"
  rm -f "$bam" "${bam}.bai"
  # Also flush GRIDSS workingdir cache to reclaim space if VCF produced
  if [[ -s "$vcf" ]]; then
    rm -rf "$dir/workingdir"
  fi
  echo "$vcf"
}

run_svaba_targeted() {
  local sample=$1 cram=$2 dir=$OUT/svaba/$sample
  local vcf=$dir/$sample.svaba.svaba.sv.vcf
  if [[ -s "$vcf" ]]; then echo "$vcf"; return 0; fi
  mkdir -p "$dir"
  # SvABA's -k wants a comma-separated region string ("chr3,chr14,..."), NOT
  # a BED file. Easiest path: re-use the pre-extracted targeted BAM Delly and
  # TIDDIT already produce — SvABA on a small BAM doesn't need -k.
  local bam
  bam=$(_prepare_targeted_bam "$sample" "$cram")
  [[ -z "$bam" || ! -s "$bam" ]] && { log "  svaba: no targeted BAM, skipping"; echo ""; return 1; }
  pushd "$dir" >/dev/null
  mamba run -n svaba svaba run \
    -t "$bam" -G "$REF" -a "$sample.svaba" -p "$THREADS" \
    >"$OUT/logs/svaba_${sample}.log" 2>&1 || true
  popd >/dev/null
  echo "$vcf"
}

run_tiddit_targeted() {
  local sample=$1 cram=$2 dir=$OUT/tiddit/$sample
  local vcf=$dir/$sample.tiddit.vcf
  if [[ -s "$vcf" ]]; then echo "$vcf"; return 0; fi
  mkdir -p "$dir"
  # Use the same pre-extracted targeted BAM as Delly (shared NAS read).
  local bam
  bam=$(_prepare_targeted_bam "$sample" "$cram")
  [[ -z "$bam" || ! -s "$bam" ]] && { log "  tiddit: no targeted BAM, skipping"; echo ""; return 1; }
  pushd "$dir" >/dev/null
  mamba run -n tiddit tiddit --sv --skip_assembly \
    --bam "$bam" --ref "$REF" -o "$sample.tiddit" \
    --threads "$THREADS" \
    >"$OUT/logs/tiddit_${sample}.log" 2>&1 || true
  [[ -f "$sample.tiddit.vcf" ]] && mv "$sample.tiddit.vcf" "$vcf" 2>/dev/null || true
  popd >/dev/null
  echo "$vcf"
}

# ----- per-sample dispatch -----------------------------------------------

IFS=',' read -ra TOOL_LIST <<< "$TOOLS"

for sample in $SAMPLES; do
  cram=${CRAM_OF[$sample]:-}
  if [[ -z "$cram" || ! -f "$cram" ]]; then
    log "SKIP $sample — no CRAM (${cram:-unset})"
    continue
  fi
  for tool in "${TOOL_LIST[@]}"; do
    log "[$tool] $sample — start"
    start=$(date +%s)
    case "$tool" in
      manta)  vcf=$(run_manta_targeted   "$sample" "$cram");;
      delly)  vcf=$(run_delly_targeted   "$sample" "$cram");;
      svaba)  vcf=$(run_svaba_targeted   "$sample" "$cram");;
      gridss) vcf=$(run_gridss_targeted  "$sample" "$cram");;
      tiddit) vcf=$(run_tiddit_targeted  "$sample" "$cram");;
      *) log "  unknown tool: $tool"; continue;;
    esac
    elapsed=$(( $(date +%s) - start ))
    log "[$tool] $sample — done in ${elapsed}s, vcf=$vcf"
  done
  # Tear down the per-sample targeted BAM after all tools for this sample
  # have run — saves /tmp space across the loop.
  _cleanup_targeted_bam "$sample"
done

log "All caller passes done. Run scripts/score_external_callers.py next."
