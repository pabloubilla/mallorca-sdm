# Mallorca SDM

Benchmark de **modelos de distribución de especies profundos** (Deep SDM) con observaciones de plantas terrestres en Mallorca, integrando varias fuentes de datos y evaluando transferencia entre tipos de muestreo.

## Fuentes y clases de observación

| Clase   | Significado                                    | Origen (aprox.)                    | Unidad en el modelo |
| ------- | ---------------------------------------------- | ---------------------------------- | ------------------- |
| **PAC** | Presencia–ausencia completa (cuadrícula Flora) | Flora + GBIF (incertidumbre 707 m) | Cuadrícula `COD10X10` |
| **PAU** | Presencia–ausencia incompleta                  | Flora sin completar cuadrícula     | Cuadrícula `COD10X10` |
| **POV** | Presencia validada (puntos GPS)                | Biodibal                           | Punto (`lon`/`lat`) |
| **POU** | Presencia no validada (puntos GPS)             | GBIF                               | Punto (`lon`/`lat`) |

Los nombres científicos se unifican contra el backbone de GBIF en `data_preparation.R`.

## Entorno Python

```bash
cd mallorca-sdm
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install torch scikit-learn matplotlib SciencePlots   # modelado, curvas, estilo nature
```

Figuras (`sampling_effort.py`, `map_grid_density_source.py`): estilo [SciencePlots nature](https://github.com/garrettj403/SciencePlots) vía `plot_style.py` (`science` + `nature` + `no-latex`).

Usa siempre el `pip` del venv (`which pip` → `.venv/bin/pip`). Dependencias base en `requirements.txt`.

## Estructura del flujo de trabajo

Ejecutar los pasos **en este orden** (cada uno depende del anterior):

```
1. data_preparation.R       → data/raw/full_data_clean_saez.gpkg
2. extract_raster_data.py   → data/raw/extracted_data.csv (+ .gpkg)
3. preprocessing.py         → data/processed/ (matrices y covariables)
4. run_model.py             → output/best_model_*.pt, output/auc_per_species_*.csv
5. sampling_effort.py       → curvas esfuerzo de muestreo (opcional)
```

### 1. `data_preparation.R` (R)

Carga Biodibal, GBIF, Flora/Sáez; filtra y fusiona geometrías; asigna `class`; resuelve taxonomía con **rgbif** (`name_backbone` / revisión manual en `taxonkey_list_review.xlsx`). Salida principal: `data/raw/full_data_clean_saez.gpkg`.

### 2. `extract_raster_data.py` (Python)

Muestrea covariables raster (`.tif`) en cada registro del GPKG. Configurar rutas en `RASTER_FOLDERS` (Bioclim, CORINE, pendiente, elevación, etc.). Deduplica geometrías para el extract y reexpande a todas las filas. Salida: `data/raw/extracted_data.csv` (y opcionalmente `.gpkg`).

### 3. `preprocessing.py` (Python)

- Join espacial a `data/raw/Flora_net.gpkg` → columna `COD10X10`.
- **PAC / PAU:** `species_matrix_{cls}.csv` — una fila por celda `UTMCODE1X1`.
- **POV / POU (punto):** mismos nombres sin sufijo — fila por GPS; bioclim en el punto.
- **POV / POU (grouped):** `species_matrix_{cls}_grouped.csv` — media bioclim y presencias por celda (como PA).
- `grid_cells_{cls}.csv` y `grid_cells_{cls}_grouped.csv` para el muestreo por celdas.

**Catálogo de especies:** presentes en ≥3 celdas PAC (`MIN_SPECIES_LIST_PAC_GRIDS`).

**Umbrales al entrenar** (`run_model.py`, sobre el dataset completo de la fuente):

| Fuente | Mínimo para incluir especie |
|--------|----------------------------|
| PAU    | ≥10 celdas con presencia   |
| POV/POU (punto) | ≥10 puntos con presencia |
| POV/POU (`_grouped`) | ≥10 celdas con presencia |

**Evaluación en PAC:** especies con ≥3 celdas PAC (catálogo).

Filas con `slope` NA se descartan. Listas de variables en `covariables.py`. Salida en `data/processed/`.

### 4. `run_model.py` (Python)

Entrena un MLP en una fuente y evalúa AUC por especie en **PAC** (cuadrículas).

| Fuente | Pérdida              | Datos de entrenamiento |
|--------|----------------------|-------------------------|
| PAU    | `BalancedBCELoss`    | Subconjunto de **celdas** |
| POV/POU | `DeepMaxEntLoss`   | Puntos dentro de celdas muestreadas, o celdas si `grouped=True` |

`size_train` = celdas muestreadas. `grouped=True` usa ficheros `_grouped` (PO alineado con PA).

```bash
python run_model.py
python -c "from run_model import run_model; run_model('POV', size_train=10)"
python -c "from run_model import run_model; run_model('POV', size_train=10, grouped=True)"
```

### 5. `sampling_effort.py` (Python, opcional)

Compara **PAU, POU y POV** en modo **point** y **grouped** (PO agregado a celda). Figura **3×3**:

| Fila | Contenido |
|------|-----------|
| 1 | AUC vs `n_grids` (modo point) |
| 2 | AUC vs `n_train_rows` (modo point) |
| 3 | AUC vs `n_grids` (modo grouped; filas ≈ celdas también en PO) |

Columnas: bins de prevalencia (% celdas PAC). CSV incluye columna `data_mode`.

**Bins de prevalencia en PAC** (para estratificar la AUC):

| Bin    | Criterio (% de celdas PAC) |
|--------|----------------------------|
| rare   | &lt; 5%                    |
| medium | 5% – 30%                   |
| common | &gt; 30%                   |

#### Conjunto de especies para la evaluación

| Modo | Descripción |
|------|-------------|
| `all` | Sin filtro extra (por defecto). |
| `intersection_pac_pov_pou` | Solo especies con ≥1 presencia en PAC, POV y POU (unidades naturales de cada matriz). |

```bash
python sampling_effort.py
python sampling_effort.py --eval-species intersection_pac_pov_pou
```

```bash
python intersection_POV_POU_PAC.py   # listar tamaño de la intersección
```

Salidas: `output/sampling_effort_auc.csv`, `output/sampling_effort_auc_by_prevalence.png` (sufijo `_intersection_pac_pov_pou` si aplica).

## Análisis y utilidades

| Archivo | Uso |
|---------|-----|
| `EDA_data.ipynb` | Exploración de datos |
| `analyse_data.R` | Mapas y resúmenes espaciales |
| `intersection_POV_POU_PAC.py` | Especies en PAC ∩ POV ∩ POU |
| `sampling_effort.py` | Curvas celdas de entrenamiento vs. AUC en PAC |
| `covariables.py` | Lista de columnas bioclim |

## Requisitos

- **R:** `sf`, `tidyverse`, `rgbif`, …
- **Python 3.12+:** `geopandas`, `rasterio`, `pandas`, `numpy`, `torch`, `scikit-learn`, `matplotlib`

Los rasters y datos crudos (`data/raw/`) no están versionados en git; `output/` está en `.gitignore`.

## Nota sobre el experimento actual

El benchmark de transferencia es **entrenar en PAU / POU / POV → evaluar en PAC** (cuadrícula Flora). El esfuerzo de muestreo se expresa en **celdas** comunes (`COD10X10`); PO conserva resolución a punto dentro de esas celdas. La intersección taxonómica (`intersection_pac_pov_pou`) es opcional para equilibrar qué especies entran en la media de AUC entre fuentes.

Con el subconjunto espacial actual del extract hay pocas decenas de celdas con datos; al ampliar cobertura en Mallorca, los porcentajes de prevalencia y los valores de `SIZE_GRID` tendrán más rango.
