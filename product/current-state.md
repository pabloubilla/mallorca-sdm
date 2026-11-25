# Current State

## Scope

Beyond PO and PA — niveles PAC/PAU/POV/POU; inventario + transferencia SDM a PAC.

## Code

GitHub: `pabloubilla/mallorca-sdm` (branch `main`).

Pipeline:

1. `data_preparation.R` / `rebuild_full_data_saez.py` → `data/raw/full_data_saez.gpkg`
2. `scripts/apply_curated_taxon_filter.py` → keep flora∪alóctonos → `full_data_saez_curated.gpkg` + `data/raw/curated/`
3. `extract_raster_data.py` → `extracted_data.csv` (filtro curated reusa extract + `pa_occurrences`)
4. `preprocessing.py --curated` → `data/processed/curated/`
5. Análisis: `scripts/analyze_rarity_by_source.py`, `scripts/species_accumulation_mix.py`
6. SDM: `sampling_effort.py` / `optimal_model.py` con `MALLORCA_PROCESSED_DIR=data/processed/curated` → `output/curated/`

## Data (curated)

- Keep: **1545** taxones (`en_clave_o_aloctonas=SI`); GPKG **1543** spp / **301 242** obs
- Matrices curated: catálogo PAC ≥1 celda → **1259** spp
- Listado alóctonos canónico: `data/raw/Alien_Flora_Balearic_Records.xlsx`
- Review: `data/raw/dataset_clave_aloctonas_review.xlsx`
- **Share pack in git** (~53 MB): `extracted_data.csv.gz` + `pa_occurrences.csv` + keep + `Flora_net.gpkg` — colaborador corre `preprocessing.py --curated` (matrices no van al repo)

## Manuscript

- Paper canónico: `manuscript/mallorca-sdm-manuscript/main-2.tex` (Beyond PO and PA)
- `main.tex` (dataset MallorcaPlant): **archivado / no primario**
- Figuras/resultados curated: `output/curated/`

## Infrastructure

- Local: `~/research/mallorca-sdm`
- NAS: `smb://dxp4800pro-14df.local/personal_folder/research/mallorca-sdm`
- Product Lead: Iván Cortés-Fernández
- Colaborador: Pablo Ubilla (INRIA Montpellier)

## Open checklist

- [x] Product Lead = Iván; Pablo = colaborador INRIA
- [x] Biel FLORA Mallorca + cuadrículas PAC/PAU
- [x] Filtro curated flora Mll ∪ alóctonos aplicado
- [ ] Validar offsets A–D 500 m con Biel
- [ ] Clusters PAC + LOCO
- [ ] Journal target
- [ ] Lock métricas finales post cobertura raster completa
- [ ] Rescates manuales 412 drop / `rescate_texto`
