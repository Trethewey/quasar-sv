#!/usr/bin/env bash
# Manual verification: how many BCL6-region reads have their mate on chr14
# (true t(3;14) signal)?
# Usage: check_bcl6_igh.sh <cram> [reference]
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/config.sh" 2>/dev/null || true
CRAM=${1:?Usage: $0 <cram> [reference]}
REF=${2:-${FF_REFERENCE:?reference required (positional or FUSIONFORGE_REFERENCE env)}}

echo "Pass 1: BCL6 reads whose mate is on chr14"
samtools view --reference "$REF" "$CRAM" chr3:187420000-187800000 2>/dev/null \
  | awk -F'\t' 'BEGIN{n=0} $7=="chr14" {n++; if (n<=8) print $3":"$4, "cigar="$6, "mate="$7":"$8} END{print "BCL6 region reads with chr14 mate:", n}'

echo
echo "Pass 2: BCL6 reads with SA tag pointing to chr14"
samtools view --reference "$REF" "$CRAM" chr3:187420000-187800000 2>/dev/null \
  | awk -F'\t' 'BEGIN{n=0} {for(i=12;i<=NF;i++){if ($i ~ /^SA:Z:chr14,/) {n++; if (n<=8) print $3":"$4, "cigar="$6, $i; break}}} END{print "BCL6 region reads with SA->chr14:", n}'

echo
echo "Pass 3: chr2:32916xxx reads whose REAL mate (mate-pos -> chr14?) "
samtools view --reference "$REF" "$CRAM" chr2:32915800-32916800 2>/dev/null \
  | awk -F'\t' 'BEGIN{n=0; total=0} {total++; if ($7=="chr14") n++} END{print "reads at chr2 artefact total:", total, "  mate=chr14:", n}'
