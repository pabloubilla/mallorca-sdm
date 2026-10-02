# Current State

## Scope

Deep SDM Mallorca — benchmark fuentes PAC/PAU/POV/POU + esfuerzo de muestreo.

## Code

GitHub: `pabloubilla/mallorca-sdm` (branch `main`).

Pipeline (orden README):

1. `data_preparation.R` / `rebuild_full_data_saez.py` → `data/raw/full_data_saez.gpkg`
   - Fuentes: GBIF Flora (~707 m), **Saez/Puig_Major** (500 m), **Biel FLORA GBO Mallorca**, Biodibal, GBIF PO
   - Biel cuadrículas → PAC vs PAU para GBO; columnas `source` + `coordinateUncertaintyInMeters` en GPKG
   - Sin covariables (extract aparte)
2. `extract_raster_data.py` → `extracted_data.csv` / `.gpkg`
3. `preprocessing.py` → matrices `data/processed/`
4. `run_model.py` / `sampling_effort.py` / `optimal_model.py` → `output/`

Diseño abierto: `notes.md` (loss mixes, hold-out PAC, CV espacial, occupancy Later).

## Data

- `data/raw/Biel/`: FLORA GBO + Quadrícules + illots
- `data/raw/Puig_Major.csv`: PA 500 m (Saez)
- `output/` corridas previas — **no métricas locked** hasta re-run post-Biel

## Manuscript

- `manuscript_draft.md` en raíz
- `manuscript/` LaTeX Wiley: **missing**

## Infrastructure

- Local: `~/research/mallorca-sdm`
- NAS: `smb://dxp4800pro-14df.local/personal_folder/research/mallorca-sdm`
- Product Lead: Iván Cortés-Fernández
- Colaborador: Pablo Ubilla (INRIA Montpellier)
- Portfolio card: `portfolio_research/projects/mallorca-sdm.md`

## Open checklist

- [x] Product Lead = Iván; Pablo = colaborador INRIA
- [x] Biel FLORA Mallorca + cuadrículas PAC/PAU en `data_preparation.R`
- [x] `source` + incertidumbre por grano en esquema merge
- [x] Re-run prep → `full_data_saez.gpkg` + taxonkey (Biel incluido; extract/preprocess pendiente)
- [ ] Validar offsets A–D 500 m con Biel
- [ ] Clusters PAC + LOCO
- [ ] Journal + LaTeX Wiley
- [ ] Lock métricas desde `output/` post re-run
