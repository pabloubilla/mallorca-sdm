# Mallorca SDM

Deep species distribution models (Deep SDM) for terrestrial plants in Mallorca.
Benchmark: train on incomplete / presence-only sources, evaluate transfer on complete presence–absence grids (**PAC**).

Collaborators who want to work on **modeling methodology** can start from the curated share pack below (no rasters required).

## Observation classes

| Class | Meaning | Typical origin | Model unit |
| ----- | ------- | -------------- | ---------- |
| **PAC** | Presence–absence, complete sampling | Flora + Sáez/Puig Major + Biel GBO (mapped cells) | Cell → centroid + uncertainty |
| **PAU** | Presence–absence, incomplete | Incomplete flora + Biel outside PAC cells | Same |
| **POV** | Validated presence (GPS) | Biodibal | Point (`lon`/`lat`) |
| **POU** | Unvalidated presence (GPS) | GBIF | Point (`lon`/`lat`) |

`PAC` / `PAU` = **sampling quality**, not cell size. Column `source` ∈ {`GBIF`, `Saez`, `Biel`, `Biodibal`}.

### Curated taxon filter

Shared data keep taxa that are in the **official Mallorca flora** *or* the **official allochthonous list** (`en_clave_o_aloctonas == SI`). Everything else is dropped as nomenclature / geography noise.

- Keep list: `data/raw/curated_taxon_keep.csv`
- Filter script: `scripts/apply_curated_taxon_filter.py`
- Review workbook (local / not in git): `data/raw/dataset_clave_aloctonas_review.xlsx`

## Quickstart (modeling on curated data)

### 1. Environment

```bash
git clone https://github.com/pabloubilla/mallorca-sdm.git
cd mallorca-sdm
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install torch scikit-learn matplotlib SciencePlots
```

### 2. Curated inputs (in git)

| Path | Role |
| ---- | ---- |
| `data/raw/curated/extracted_data.csv.gz` | Raster pixels / GPS points + covariates |
| `data/raw/curated/pa_occurrences.csv` | PA labels (species × 1 km cell) |
| `data/raw/curated/filter_log.txt` | Filter counts |
| `data/raw/curated_taxon_keep.csv` | Taxa kept |
| `data/raw/Flora_net.gpkg` | 1 km grid (`COD10X10` for PO) |

Species **matrices are not shipped**. Build them locally:

```bash
python preprocessing.py --curated
# → data/processed/curated/{species_matrix,covariates,grid_cells}_*.csv
```

`preprocessing.py` accepts `.csv` or `.csv.gz` for the extracted table.

### 3. Train / evaluate

```bash
export MALLORCA_PROCESSED_DIR=data/processed/curated

python -c "from run_model import run_model; run_model('PAU', size_train=10)"
python -c "from run_model import run_model; run_model('POV', size_train=10, grouped=True)"

python sampling_effort.py
python optimal_model.py
```

Canonical fair mixes use **cell grain** (`grouped=True` for PO; PA pixels aggregated to `UTMCODE1X1` for AUC). Outputs default under `output/` (gitignored); curated analysis scripts often write `output/curated/`.

| Source | Loss | Training rows |
| ------ | ---- | ------------- |
| PAU | `BalancedBCELoss` | Sampled cells (pixel rows inside cells) |
| POV / POU | `DeepMaxEntLoss` | Points in sampled cells, or cells if `grouped=True` |

Evaluation: AUC per species on **PAC**.

## Full pipeline (optional — needs local raw + rasters)

Raw GBIF/Biodibal/Biel stacks and rasters are **not** in git (`data/` ignored except the curated pack above).

```
1. data_preparation.R / rebuild_full_data_saez.py  → data/raw/full_data_saez.gpkg
2. scripts/apply_curated_taxon_filter.py           → data/raw/curated/ + *_curated.gpkg
3. extract_raster_data.py                         → extracted_data (+ curated reuse)
4. preprocessing.py [--curated]                   → data/processed/[curated/]
5. run_model.py / sampling_effort.py / optimal_model.py
```

### Preprocessing notes

- Join to `Flora_net.gpkg` → `COD10X10` on PO.
- **PAC / PAU:** one row per ~100 m raster pixel; species labels at cell `UTMCODE1X1`.
- **POV / POU:** point matrices; optional `*_grouped.csv` (cell mean).
- Species catalog: PAC ≥ `MIN_SPECIES_LIST_PAC_GRIDS` cells (curated default = 1).
- Drop rows with covariate NA (e.g. `slope`). Columns in `covariables.py`.

### Sampling effort (optional)

```bash
python sampling_effort.py
python sampling_effort.py --eval-species intersection_pac_pov_pou
```

Prevalence bins on PAC: rare &lt;5%, medium 5–30%, common &gt;30%.

## Requirements

- **Python 3.12+:** see `requirements.txt`; plus `torch`, `scikit-learn`, `matplotlib` for modeling
- **R** (full prep only): `sf`, `tidyverse`, `rgbif`, …

## Main scripts (light map)

### Pipeline / modeling (repo root)

| Script | Role |
| ------ | ---- |
| `data_preparation.R` / `rebuild_full_data_saez.py` | Build `full_data_saez.gpkg` from raw GBIF / Biel / Biodibal / Sáez |
| `extract_raster_data.py` | Sample rasters → `extracted_data` (+ PA pixel rows) |
| `preprocessing.py` | Build `species_matrix_*`, `covariates_*`, `grid_cells_*` |
| `covariables.py` | Covariate column list |
| `models.py` | MLP + `BalancedBCELoss` / `DeepMaxEntLoss` (+ phylo Laplace helper) |
| `run_model.py` | Train one source → per-species AUC on PAC |
| `sampling_effort.py` | Learning curves vs effort / prevalence bins |
| `optimal_model.py` | Factorial PAU/POV/POU mixes → mean PAC AUC |
| `plot_style.py` | SciencePlots nature style for figures |
| `map_grid_density_source.py` | 1 km density maps by observation class |

### Analysis / prep helpers (`scripts/`)

| Script | Role |
| ------ | ---- |
| `apply_curated_taxon_filter.py` | Keep flora ∪ allochthonous → curated GPKG + extract pack |
| `build_dataset_clave_aloctonas_review.py` | Review workbook clave ↔ dataset ↔ allochthonous |
| `extract_clave_vs_dataset.py` | Parse clave / match vs dataset names |
| `analyze_rarity_by_source.py` | Rarity bins + PO “rarification” vs PAC |
| `species_accumulation_mix.py` | Inventory accumulation under source mixes |
| `analyze_taxonomy_by_source.py` | Family/order composition PAC/PAU/POV/POU |
| `build_opentree_laplacian.py` | OpenTree TNRS + cophenetic Laplacian cache |
| `run_phylo_laplace_ab.py` | Smoke: MLP ± phylo Laplace → PAC AUC |
| `plot_illustrative_beyond_po_pa.py` | Illustrative manuscript figures |

Product docs: `product/` (vision, decisions, roadmap, current-state). Manuscript: `manuscript/mallorca-sdm-manuscript/` (canonical `main-2.tex`).

## Experiment note

Transfer benchmark: train **PAU / POU / POV** → evaluate **PAC**. Sampling effort is measured in shared cells; PO can stay at point resolution inside those cells unless `grouped=True`.
