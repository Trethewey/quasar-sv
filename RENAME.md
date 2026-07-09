# Renaming the package

> `quasarsv` is a **working name**. Final name is undecided.

When you settle on the final name, do the rename in one command:

```bash
# Dry-run first to see what will change
bash scripts/rename_package.sh <new_name> --dry-run

# Apply it
bash scripts/rename_package.sh <new_name>
```

## What the rename script does

1. Renames `src/quasarsv/` → `src/<new_name>/`
2. Rewrites every case-sensitive occurrence of `quasarsv` in code, docs,
   scripts, tests, `pyproject.toml`, `LICENSE`, `README.md`, `RESUME_NOTE.md`
3. Re-runs the test suite (currently 64 tests) — aborts if anything breaks

## What the script does NOT touch (do these by hand if needed)

* `archive/` — historic snapshots; the old name is a feature there
* `output/` — fully regeneratable from source, no need to touch

## Constraints on the new name

The script enforces:

* lowercase letters, digits, underscores only
* starts with a letter (Python package identifier rules)
* not a Python reserved keyword

It does NOT check:

* PyPI availability (do `pip search <name>` or visit pypi.org/project/<name>/)
* GitHub repo availability
* Trademark conflicts (be careful — `fusioncatcher`, `STAR-Fusion`, `Arriba`,
  `FusionInspector`, `Manta`, `Delly`, `SvABA`, `GRIDSS`, `TIDDIT` are all
  taken; pick something distinct)

## Files that intentionally reference the working name

These deliberately preserve the `quasarsv` text and SHOULD NOT change
when you rename:

* `archive/session_logs/*` — historical snapshots
* `archive/competitor_tools/factera_distribution/` — competitor tool, unrelated

After the rename, also regenerate the DOCX vignette so its embedded
content uses the new name:

```bash
python3 scripts/make_vignette_docx.py
```
