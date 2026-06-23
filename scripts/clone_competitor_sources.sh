#!/usr/bin/env bash
# Clone source repos of each competitor SV caller for source-code survey.
set -euo pipefail
DST=${1:-$HOME/tools_src}
mkdir -p "$DST"
cd "$DST"
for repo in \
  Illumina/manta \
  dellytools/delly \
  walaj/svaba \
  arq5x/lumpy-sv \
  SciLifeLab/TIDDIT \
  PapenfussLab/gridss \
  AstraZeneca-CGI/svict \
  stanford-sentinel/factera; do
  name=${repo##*/}
  if [ -d "$name" ]; then
    echo "[skip] $name (already cloned)"
    continue
  fi
  echo "[clone] $repo"
  git clone --depth 1 "https://github.com/$repo.git" "$name" 2>&1 | tail -2 || \
    echo "  -> clone failed for $repo"
done
ls -la "$DST"
