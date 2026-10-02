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

## Open (USER)

- SoT datos: clone local vs NAS DX
- Journal target
- Convenio exacto A/B/C/D de offsets 500 m (ahora A=SW, B=SE, C=NW, D=NE) — validar con Biel
- Umbrales de inclusión de especies (PAC ≥3 celdas, etc.)
- Diseño final CV espacial (n clusters, LOCO vs k-fold)
- Data Availability (GBIF / Flora / Biodibal / GBO)
- Reposicionar centroides en celdas muy marinas / illots → Discussion Later (no en prep Now)
- Criterio “gran part mar” → Discussion Later (no define PAC en código)
