# Decisions

## Locked

- Área: **Mallorca only** (Menorca / Eivissa fuera por ahora)
- Taxón: plantas terrestres vasculares
- Product Lead: **Iván Cortés-Fernández**; colaborador: Pablo Ubilla Pavez (INRIA Montpellier / EVERGREEN)
- Cuatro clases: **PAC**, **PAU**, **POV**, **POU** (calidad de muestreo; **PAC ≠ resolución**)
- PA: centroide de celda + `coordinateUncertaintyInMeters` = semidiagonal (1 km ≈ 707 m; 500 m ≈ 354 m)
- Fuentes en GPKG: `source` ∈ {GBIF, Saez/Puig_Major, Biel, Biodibal}
- Biel FLORA GBO Mallorca entra en prep; **cuadrículas GBO** definen qué celdas Biel son PAC vs PAU (mapa flora)
- Puig Major (`Puig_Major.csv`): PA 500 m (Saez), unc ≈ 354 m
- Evaluación canónica: **PAC** (AUC por especie; media)
- PA extract: un píxel raster (~100 m) por fila dentro de celda 1×1 km; agregar predicciones a celda para AUC
- PO: GPS punto; opcional `_grouped` a celda para comparar con PA
- Loss: PA → `BalancedBCELoss`; solo PO → `DeepMaxEntLoss`
- Stack: R (`data_preparation.R`) + Python 3.12 + PyTorch
- Remote: `https://github.com/pabloubilla/mallorca-sdm`
- Integrity: no inventar números; no autores/títulos en artwork
- **Curated taxon set (2026-11-16):** excluir del dataset curado los taxones que **no** estén en la flora de Mallorca (clave / Quarta edició, `clave_Mll`) **ni** en el listado oficial de taxones alóctonos. Interpretación: errores / ruido de nomenclatura o geografía → mejor no incluir. Keep si está en flora Mallorca **o** en listado alóctonos.
- **Listado alóctonos canónico:** `data/raw/Alien_Flora_Balearic_Records.xlsx` (join en `scripts/build_dataset_clave_aloctonas_review.py`; filtro en `scripts/apply_curated_taxon_filter.py`).
- **Paper framing (2026-11-25):** canónico = *Beyond PO and PA* (niveles PAC/PAU/POV/POU); paper-dataset GeoPlant/MallorcaPlant descartado como primario.
- **Mixes SDM justos:** PO y PA a grano de celda (`grouped=True`; PA pixels agregados a `UTMCODE1X1`).
- Trabajo fino clave↔GBIF (`gbif_status`, reconciliación nombre a nombre en `clave_mallorca_especies.xlsx` / REVISION): **stand-by**; no bloquea el filtro de exclusión anterior.
- Rescates `rescate_texto` / revisión 412 drop: **Later** (no bloquean curated Now).

## Open (USER)

- SoT datos: clone local vs NAS DX
- Journal target
- Convenio exacto A/B/C/D de offsets 500 m (ahora A=SW, B=SE, C=NW, D=NE) — validar con Biel
- Umbrales de inclusión de especies (PAC ≥3 celdas, etc.)
- Diseño final CV espacial (n clusters, LOCO vs k-fold)
- Data Availability (GBIF / Flora / Biodibal / GBO)
- Reposicionar centroides en celdas muy marinas / illots → Discussion Later (no en prep Now)
- Criterio “gran part mar” → Discussion Later (no define PAC en código)
