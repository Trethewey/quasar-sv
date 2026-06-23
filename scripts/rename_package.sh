#!/usr/bin/env bash
# rename_package.sh — change the working name "quasarsv" to anything else.
#
# What it does:
#   1. Renames src/quasarsv/   ->  src/<newname>/
#   2. Rewrites every "quasarsv" occurrence in code, scripts, docs, tests,
#      pyproject.toml, LICENSE, README, RESUME_NOTE — case-sensitive only.
#   3. Reruns the test suite to confirm nothing broke.
#
# What it does NOT touch:
#   * Files inside `archive/`   (historic; deliberate snapshot)
#   * Files inside `output/`    (regeneratable artefacts)
#   * Memory entries (~/.claude/projects/.../memory/) — those describe
#     facts about *this* project and you should update them by hand if the
#     name change is permanent
#
# Usage:   bash scripts/rename_package.sh <new_name>
# Example: bash scripts/rename_package.sh ivyforge
#
# IMPORTANT — the new name must be a valid Python package identifier:
#   * lowercase letters, digits, underscores only
#   * starts with a letter
#   * not a Python keyword
#
# Dry-run first: bash scripts/rename_package.sh <new_name> --dry-run

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OLD_NAME="quasarsv"

NEW_NAME="${1:-}"
DRY_RUN=0
[[ "${2:-}" == "--dry-run" ]] && DRY_RUN=1

if [[ -z "$NEW_NAME" ]]; then
  echo "Usage: $0 <new_name> [--dry-run]"
  exit 2
fi

if ! [[ "$NEW_NAME" =~ ^[a-z][a-z0-9_]*$ ]]; then
  echo "ERROR: '$NEW_NAME' is not a valid Python package identifier" >&2
  echo "  Required: lowercase letters/digits/underscores, starting with a letter" >&2
  exit 2
fi

PYTHON_KEYWORDS="False None True and as assert async await break class continue def del elif else except finally for from global if import in is lambda nonlocal not or pass raise return try while with yield"
for kw in $PYTHON_KEYWORDS; do
  [[ "$NEW_NAME" == "$kw" ]] && { echo "ERROR: '$NEW_NAME' is a Python keyword" >&2; exit 2; }
done

echo "Rename plan:"
echo "  $OLD_NAME -> $NEW_NAME"
echo "  src/$OLD_NAME/ -> src/$NEW_NAME/"
echo "  pyproject.toml: name + console-script entry"
echo "  Code, docs, scripts, tests: every occurrence (case-sensitive)"
echo

if [[ $DRY_RUN -eq 1 ]]; then
  echo "[dry-run] Files that would be modified:"
  grep -rl --include='*.py' --include='*.sh' --include='*.md' --include='*.toml' \
    --include='*.tsv' --include='*.txt' \
    --exclude-dir=archive --exclude-dir=output --exclude-dir='__pycache__' \
    --exclude-dir='.pytest_cache' \
    "$OLD_NAME" "$PROJECT_ROOT"
  echo
  echo "[dry-run] Done. Re-run without --dry-run to apply."
  exit 0
fi

# Safety: confirm
read -p "Type the new name '$NEW_NAME' to confirm rename: " CONFIRM
if [[ "$CONFIRM" != "$NEW_NAME" ]]; then
  echo "ABORTED — confirmation mismatch."
  exit 3
fi

echo "[1/3] Renaming src/$OLD_NAME/ -> src/$NEW_NAME/"
if [[ -d "$PROJECT_ROOT/src/$OLD_NAME" ]]; then
  mv "$PROJECT_ROOT/src/$OLD_NAME" "$PROJECT_ROOT/src/$NEW_NAME"
fi

echo "[2/3] Rewriting string references in code/docs/scripts/tests/configs"
TARGETS=$(grep -rl --include='*.py' --include='*.sh' --include='*.md' \
  --include='*.toml' --include='*.tsv' --include='*.txt' --include='LICENSE' \
  --exclude-dir=archive --exclude-dir=output --exclude-dir='__pycache__' \
  --exclude-dir='.pytest_cache' --exclude-dir='.claude' \
  "$OLD_NAME" "$PROJECT_ROOT" 2>/dev/null || true)
if [[ -n "$TARGETS" ]]; then
  echo "$TARGETS" | while read -r f; do
    [[ -f "$f" ]] && sed -i "s|$OLD_NAME|$NEW_NAME|g" "$f"
  done
  echo "    $(echo "$TARGETS" | wc -l) files updated"
fi

echo "[3/3] Re-running test suite to confirm rename is clean"
cd "$PROJECT_ROOT"
PYTHONPATH=src python3 -m pytest tests/ -q || {
  echo "TESTS FAILED — investigate before committing." >&2
  exit 4
}

echo
echo "Rename complete: $OLD_NAME -> $NEW_NAME"
echo "Next steps:"
echo "  - Manually review: README.md, RESUME_NOTE.md, docs/*.md"
echo "  - Re-render the docx vignette:  python3 scripts/make_vignette_docx.py"
echo "  - Update auto-memory entries if the name change is permanent"
echo "  - Recreate the bg job state if any bg shells still reference '$OLD_NAME'"
