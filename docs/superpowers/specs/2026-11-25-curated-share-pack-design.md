# Curated share pack for modeling collaborators

Date: 2026-11-25  
Status: approved

## Goal

Let a collaborator clone the repo, obtain the post-filter curated occurrence tables, build matrices locally, and work on SDM methodology — without rasters or the full raw stack.

## Keep rule

Taxa with `en_clave_o_aloctonas == SI` in `dataset_clave_aloctonas_review.xlsx` (flora oficial Mallorca ∪ listado alóctonos). Materialized by `scripts/apply_curated_taxon_filter.py`.

## In-repo data (git exceptions under `data/`)

| Path | Role |
|------|------|
| `data/raw/curated/extracted_data.csv.gz` | Pixels/points + covariates (~44 MB) |
| `data/raw/curated/pa_occurrences.csv` | PA species × cell labels |
| `data/raw/curated/filter_log.txt` | Filter counts |
| `data/raw/curated_taxon_keep.csv` | Keep list |
| `data/raw/Flora_net.gpkg` | 1 km grid for PO `COD10X10` |

Matrices under `data/processed/` stay gitignored; collaborator runs `preprocessing.py --curated`.

## Code

- `preprocessing.py`: resolve `extracted_data.csv` or `.csv.gz`; pandas reads gzip.
- `.gitignore`: ignore `data/*` with negation for the files above.
- README: collaborator quickstart first; full raw pipeline secondary.

## Out of scope

Push, Git LFS, shipping species matrices, rasters, full GPKG.
