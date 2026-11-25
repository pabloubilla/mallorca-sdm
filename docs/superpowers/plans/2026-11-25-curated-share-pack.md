# Plan: curated share pack

## Done

1. Spec: `docs/superpowers/specs/2026-11-25-curated-share-pack-design.md`
2. `extracted_data.csv.gz` (~44 MB)
3. `.gitignore` exceptions for share pack + `Flora_net.gpkg`
4. `preprocessing.resolve_table_path` + `--curated`
5. `apply_curated_taxon_filter.py` also writes `.csv.gz`
6. README collaborator quickstart
7. `product/current-state.md` note

## Verify

- `python preprocessing.py --curated` with only `.gz` present (or both)
- Files &lt; 100 MB each; pack ~53 MB
- `git add` share paths not ignored

## Not done (user must ask)

- Commit / push
