#!/usr/bin/env bash
# Print sample CRAM sizes from the WGS cohort.
set -u
for s in SRR1236466 SRR1236474 SRR1236478 ERR9128954 ERR9188549; do
  f=$(find "${FF_CRAM_ROOT:?set FUSIONFORGE_CRAM_ROOT}" -name "${s}.cram" 2>/dev/null | head -1)
  if [ -n "$f" ]; then
    sz_bytes=$(stat -c '%s' "$f" 2>/dev/null)
    sz_gb=$(awk -v b="$sz_bytes" 'BEGIN { printf "%.1f", b/1024/1024/1024 }')
    echo "$s ${sz_gb}GB $f"
  fi
done
