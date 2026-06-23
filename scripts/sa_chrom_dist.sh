#!/usr/bin/env bash
# For reads at the chr2:32916xxx artefact locus, print the distribution of
# their SA-tag chromosomes — these are the TRUE partners that the locus
# is absorbing.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh" 2>/dev/null || true
CRAM=${1:?Usage: $0 <cram> [reference]}
REF=${2:-${FF_REFERENCE:?reference required (positional or FUSIONFORGE_REFERENCE env)}}
echo "SA-tag chromosome distribution for reads at chr2:32915800-32916800 in $CRAM"
samtools view --reference "$REF" "$CRAM" chr2:32915800-32916800 2>/dev/null \
  | head -100000 \
  | awk -F'\t' '{
      for (i=12; i<=NF; i++) {
        if (substr($i,1,5) == "SA:Z:") {
          # strip "SA:Z:" and take first comma-separated chrom
          sa = substr($i, 6)
          split(sa, parts, ",")
          print parts[1]
          next
        }
      }
    }' \
  | sort | uniq -c | sort -rn | head -15
