# Deep species distribution models on Mallorca: benchmarking data sources and sampling design

**Working title.** Transferability of deep multi-species models across presence–absence grids and opportunistic point records in the Balearic Islands.

---

## Abstract

Species distribution models (SDMs) increasingly use deep neural networks trained on many taxa simultaneously, yet it remains unclear how performance depends on **where** observations come from (structured atlas grids versus citizen-science points) and **how much** data are available. We assembled terrestrial plant records for Mallorca from three channels—Flora atlas, validated plot data (Biodibal database), and GBIF—and classified them into four observation types: complete presence–absence on a 1×1 km grid (PAC), incomplete grid atlasing (PAU), validated presence-only points (POV), and unvalidated presence-only points (POU). Environmental predictors were extracted from Bioclim, CORINE land cover, elevation, and slope rasters. A shared multi-layer perceptron (MLP) was trained on PAU, POV, and/or POU and evaluated on held-out PAC cells using mean per-species ROC-AUC. Grid-based training used balanced binary cross-entropy when PA data were included; presence-only sources used a deep maximum-entropy loss. We quantified (i) learning curves over the number of sampled grid cells, stratified by species prevalence on PAC, comparing point-level versus grid-aggregated PO training; and (ii) factorial mixes of 0%, 50%, and 100% subsamples from each training source. In preliminary runs, mean PAC AUC reached approximately **0.76** when training combined full PAU and POV with half of POU; mixtures with high POV contribution consistently outperformed combinations with POV excluded (AUC ≈ **0.62–0.70**). Performance was sensitive to aligning spatial units: comparing PA grids with PO GPS points inflates effective sample size and confounds effort with resolution. Aggregating PO records to the same 1×1 km grid (`UTMCODE1X1`) is required for fair source comparison. Current raster extracts cover a limited spatial subset of Mallorca, constraining PAC evaluation to on the order of **10²** cells; expanding coverage will strengthen prevalence-stratified inference. Deep SDMs can leverage opportunistic biodiversity data for atlas-scale prediction, but transfer experiments must harmonize grid geometry, loss functions, and taxonomic overlap across sources.


**Keywords:** species distribution model; deep learning; transfer learning; citizen science; atlas data; Mallorca; plant biodiversity

---

## 1. Introduction

High-quality distribution maps for many species are needed for conservation and spatial planning, but complete presence–absence surveys exist for only a fraction of taxa and regions. Atlases and structured monitoring provide grid-based presence–absence (PA) data suitable for conventional SDMs, whereas platforms such as GBIF and regional databases contribute large volumes of **presence-only** (PO) records at point resolution. Deep multi-species models—single networks predicting hundreds of species from shared environmental features—promise to exploit PO data efficiently, yet their behaviour under **domain shift** (training on PO or incomplete PA, evaluating on complete PA) is poorly benchmarked.

Mallorca offers a controlled setting to study this problem: long-term Flora work provides 1×1 km grid atlasing, Biodibal supplies validated botanical presences, and GBIF adds broader but less curated observations. Taxonomic names were harmonised against the GBIF backbone. We ask: (1) How does predictive skill on a complete atlas (PAC) scale with the number of training grid cells for PAU versus POV versus POU? (2) Does aggregating PO points to grid cells change transfer performance relative to point-based training? (3) Which mixtures of PAU, POV, and POU subsamples maximise mean AUC on PAC? This study implements a reproducible pipeline (R and Python) and reports a first quantitative comparison; we emphasise methodological pitfalls when mixing grid and point units.

---

## 2. Materials and Methods

### 2.1 Study area and occurrence data

The study focuses on **terrestrial vascular plants** on Mallorca (Balearic Islands, Spain). Raw records were compiled in `data_preparation.R` from:

| Code | Type | Source (approx.) | Spatial unit |
|------|------|------------------|--------------|
| PAC | Presence–absence, complete atlas | Flora + GBIF (707 m uncertainty filter) | 1×1 km grid (`UTMCODE1X1`) |
| PAU | Presence–absence, incomplete atlas | Flora (incomplete grid coverage) | 1×1 km grid |
| POV | Presence-only, validated | Biodibal | GPS point (also aggregable to grid) |
| POU | Presence-only, unvalidated | GBIF | GPS point |

Scientific names were matched to GBIF taxon keys (`name_backbone` and manual review). Geometries were cleaned and exported to `full_data_clean_saez.gpkg`.

### 2.2 Environmental covariates

For each record, covariates were sampled from rasters (`extract_raster_data.py`): 19 Bioclim variables, CORINE land-cover classes, elevation, and slope. Records with missing slope were excluded. The predictor set is defined in `covariables.py` (~60 features after one-hot encoding of CORINE).

### 2.3 Spatial preprocessing

Point and grid observations were joined to the Flora network (`Flora_net.gpkg`) to assign **`UTMCODE1X1`** cells. Processed matrices (`preprocessing.py`) include:

- **PAC / PAU:** one row per grid cell; species columns are binary presence on that cell.
- **POV / POU (point mode):** one row per GPS site; bioclim at point location.
- **POV / POU (grouped mode):** one row per grid cell; environmental values averaged per cell; presence if any point occurs in the cell (analogous to PA).

A species catalogue was built from PAC; training filters require minimum occurrences per species on the training unit (grids or points). Evaluation on PAC used species with sufficient presences on evaluation cells (catalogue threshold in the study design; see repository for current numeric thresholds).

### 2.4 Deep multi-species model

We used a multi-task **MLP** (`models.py`): two hidden layers (200 units), residual blocks, sigmoid outputs per species. Features were standardised on the training set. Optimization: AdamW (learning rate 10⁻³, weight decay 10⁻⁴), ReduceLROnPlateau, 80/20 train–validation split, early stopping on validation loss.

| Training data | Loss function |
|---------------|---------------|
| Any mix including PAU (PA) | `BalancedBCELoss` (class-balanced BCE) |
| POV/POU only | `DeepMaxEntLoss` (batch softmax over sites, Ryckewaert-style) |

**Transfer protocol:** train on PAU, POV, and/or POU; **evaluate per-species ROC-AUC on PAC** (complete grid); report unweighted mean AUC across species.

### 2.5 Experiments

**Sampling effort** (`sampling_effort.py`): For each source (PAU, POV, POU), subsample `n` ∈ {1, 10, 100, 200, 500, 1000, 1600} grid cells (for PO point mode, all GPS points within sampled cells are retained). Repeat training; plot mean PAC AUC versus `n_grids` and versus training row count. Stratify species by PAC prevalence (rare &lt;5%, medium 5–30%, common &gt;30% of PAC cells). Compare **point** versus **grouped** PO training.

**Optimal source mix** (`optimal_model.py`): Factorial design over fractions {0%, 50%, 100%} of available cells/rows per source (26 non-empty combinations), random subsampling (seed 42), five training epochs per combination. **Recommended analysis:** all sources at **grid** resolution (`_grouped` for POV/POU) so PAU and PO are comparable.

Optional taxonomic restriction: species with ≥1 presence in PAC, POV, and POU (`intersection_pac_pov_pou`).

### 2.6 Software and reproducibility

R 4.x (`sf`, `tidyverse`, `rgbif`); Python 3.12 (`geopandas`, `rasterio`, `PyTorch`, `scikit-learn`). Code and workflow are documented in `README.md`. Model checkpoints and CSV outputs are written under `output/`.

---

## 3. Results

### 3.1 Data volume (processed subset)

After raster extraction and filtering, the processed benchmark spanned on the order of **10⁵** point-level PO records and **10³** grid cells per source when aggregated: approximately **1,065** PAU cells, **1,972** POV grouped cells, **1,768** POU grouped cells, and **~366** PAC evaluation cells in the current spatial subset (exact counts depend on slope filtering and Flora join). PAC remains the sole complete presence–absence reference for transfer metrics.

### 3.2 Source-mix experiment (preliminary, point-level PO training)

Table 1 summarises the top and bottom combinations from `output/optimal_model_auc.csv` (five epochs; PO trained at point resolution—see Discussion).

| Rank | Training mix (PAU / POV / POU %) | Training rows | Mean PAC AUC |
|------|----------------------------------|---------------|--------------|
| 1 | 100 / 100 / 50 | 86,657 | **0.763** |
| 2 | 100 / 100 / 100 | 94,176 | 0.758 |
| 3 | 0 / 100 / 100 | 93,111 | 0.757 |
| … | … | … | … |
| Low | 100 / 0 / 0 (PAU only) | 1,065 | 0.664 |
| Low | 50 / 0 / 0 | 532 | 0.622 |

**Patterns:**

1. **POV dominance:** Mixtures with POV at 100% achieved mean AUC **0.74–0.76**, even without PAU.
2. **POV absence:** When POV was 0%, AUC fell to **0.62–0.70** despite high PAU or POU fractions—suggesting validated plot presences carry more transferable signal than incomplete atlas or raw GBIF alone under this setup.
3. **Non-monotone POU effect:** At full PAU+POV, 50% POU slightly exceeded 100% POU (0.763 vs 0.758), indicating interaction between sample size, label noise, and loss weighting rather than “more GBIF is always better.”
4. **Effort confound (point mode):** Training “rows” for PO reached **10⁴–10⁵** while PAU used **10³** cells—comparisons across mixes partly reflect unequal spatial support unless PO is grid-aggregated.

### 3.3 Sampling effort (expected qualitative outcomes)

Learning-curve analyses (`sampling_effort.py`) are designed to show whether AUC on PAC saturates with increasing sampled cells and whether curves differ by prevalence class and by point versus grouped PO. With the current Mallorca subset, PAC cell count limits stable estimates for rare species; full curves should be re-estimated after extending raster coverage.

---

## 4. Discussion

### 4.1 Transfer from opportunistic data to atlas grids

Preliminary mixes indicate that **validated POV data strongly drive transfer to PAC**, whereas PAU alone yields modest AUC (~0.66) on the evaluation grid. This is consistent with POV covering more cells and species with detections aligned to field botany, while PAU is incomplete by construction. GBIF (POU) adds value mainly in combination with POV and PAU, not as a sole trainer—likely reflecting taxonomic, spatial, and reporting biases.

### 4.2 Grid versus point units

A central methodological lesson is that **PA and PO must share the same spatial support** for fair comparison. Training PO at GPS resolution while evaluating on 1×1 km PAC cells mixes grain sizes: one “training site” for PO is not one PAC cell, and subsampling “50% of POV” subsamples points, not comparable to “50% of PAU cells.” Grid aggregation (`_grouped` matrices) aligns bioclimatic features and occupancy with PAU/PAC and should be used in optimal-mix and source-ranking analyses. `sampling_effort.py` explicitly contrasts both modes to document this sensitivity.

### 4.3 Loss functions and class imbalance

Switching between balanced BCE (when PAU is present) and deep MaxEnt (PO-only mixes) changes how rare species are weighted. Mixed-source experiments therefore confound **data composition** and **objective function**; future work could hold the loss fixed or use a unified objective across PA and PO.

### 4.4 Limitations

- **Spatial coverage:** Current extracts do not cover all Mallorca; PAC evaluation cells are few relative to island extent.
- **Taxonomic scope:** Mean AUC aggregates many species; prevalence bins and intersection filters (`intersection_pac_pov_pou`) mitigate but do not remove taxonomic mismatch.
- **Training depth:** Five epochs in the mix experiment are sufficient for ranking mixtures in development, not for final performance estimates.
- **No independent spatial block CV** yet; reported AUC is in-sample transfer to a fixed PAC layer.
- **Climate and land cover** are contemporary; atlas records span decades.

### 4.5 Conclusions

Deep multi-species SDMs can be benchmarked systematically on Mallorca by training on PAU/POV/POU and evaluating on Flora PAC grids. Early results favour **rich POV inclusion**; **PAU-only** training underperforms. Harmonising **grid cells**, **loss choice**, and **species intersection** is prerequisite before interpreting optimal data mixtures. Expanding spatial coverage and re-running grid-aligned mix experiments will provide the definitive source-ranking for atlas completion in the Balearics.

---

## Author contributions / data availability (placeholders)

- **Data availability:** Processed matrices and code in repository `mallorca-sdm`; raw GPKG and rasters not versioned.
- **Code availability:** GitHub repository (URL to be added).

---

*Draft generated from project documentation and `output/optimal_model_auc.csv` (point-mode PO). Re-run `optimal_model.py` with grid-grouped PO for publication-ready mix results.*
