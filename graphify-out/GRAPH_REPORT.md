# Graph Report - mallorca-sdm  (2026-09-25)

## Corpus Check
- 22 files · ~13,605 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 238 nodes · 403 edges · 16 communities
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS · INFERRED: 2 edges (avg confidence: 0.9)
- Token cost: 0 input · 0 output

## Graph Freshness
- Built from commit: `dfe24dab`
- Run `git rev-parse HEAD` and compare to check if the graph is stale.
- Run `graphify update .` after code changes (no API cost).

## Community Hubs (Navigation)
- preprocessing.py
- run_model.py
- optimal_model.py
- extract_raster_data.py
- sampling_effort.py
- models.py
- 2. Materials and Methods
- Mallorca SDM
- Notas de diseño experimental
- Current State
- Market stress / idoneidad — mallorca-sdm
- mallorca-sdm
- Vision
- Roadmap
- Decisions

## God Nodes (most connected - your core abstractions)
1. `run_model()` - 18 edges
2. `_load_source_subset()` - 17 edges
3. `MLP` - 10 edges
4. `main()` - 10 edges
5. `grid_cells_path()` - 10 edges
6. `_load_pac_eval()` - 10 edges
7. `main()` - 9 edges
8. `_mean_auc_on_pac()` - 9 edges
9. `species_matrix_path()` - 9 edges
10. `_save_pa_pixel_products()` - 9 edges

## Surprising Connections (you probably didn't know these)
- `_run_training_epoch()` --uses--> `MLP`  [INFERRED]
  run_model.py → models.py
- `main()` --calls--> `assign_cod10_grid()`  [EXTRACTED]
  map_grid_density_source.py → preprocessing.py
- `_mean_auc_on_pac()` --calls--> `BalancedBCELoss`  [EXTRACTED]
  optimal_model.py → models.py
- `run_model()` --calls--> `BalancedBCELoss`  [EXTRACTED]
  run_model.py → models.py
- `_mean_auc_on_pac()` --calls--> `DeepMaxEntLoss`  [EXTRACTED]
  optimal_model.py → models.py

## Import Cycles
- None detected.

## Communities (16 total, 0 thin omitted)

### Community 0 - "preprocessing.py"
Cohesion: 0.25
Nodes (18): assign_cod10_grid(), _broadcast_cell_labels_to_pixels(), _cell_species_matrix(), covariates_path(), file_suffix(), grid_cells_path(), main(), _point_site_id() (+10 more)

### Community 1 - "run_model.py"
Cohesion: 0.18
Nodes (24): DataLoader, Optimizer, is_po_source(), _aggregate_pixel_scores_to_cells(), _filter_species_by_grid_occurrence(), _filter_species_by_occurrence(), _is_pa_pixel_matrix(), _load_pac_eval() (+16 more)

### Community 2 - "optimal_model.py"
Cohesion: 0.14
Nodes (23): main(), _mallorca_limits(), GeoDataFrame, Index, 2×2 map: observation density per UTMCODE1X1 cell for PAC, PAU, POV, POU., Extent of grid cells with data, padded (~Mallorca focus)., _baseline_pcts(), _build_training_mix() (+15 more)

### Community 3 - "extract_raster_data.py"
Cohesion: 0.22
Nodes (21): _assign_pa_grids(), _cell_class_map(), collect_raster_files(), _count_inside_bounds(), extract_pa_pixels(), extract_po_points(), main(), _maybe_subsample() (+13 more)

### Community 4 - "sampling_effort.py"
Cohesion: 0.20
Nodes (13): intersection_pac_pov_pou(), Index, Utilidad: especies con presencia en PAC, POV y POU (misma regla que…, Namespace, _output_paths(), _parse_args(), plot_sampling_effort(), DataFrame (+5 more)

### Community 5 - "models.py"
Cohesion: 0.06
Nodes (21): ABNLoss, ABNModel, BalancedBCELoss, BernoulliFromLogRateLoss, BiasL2Penalty, DeepMaxEntLoss, DeepMaxentLossBias, IntegratedLoss (+13 more)

### Community 6 - "2. Materials and Methods"
Cohesion: 0.09
Nodes (21): 1. Introduction, 2.1 Study area and occurrence data, 2.2 Environmental covariates, 2.3 Spatial preprocessing, 2.4 Deep multi-species model, 2.5 Experiments, 2.6 Software and reproducibility, 2. Materials and Methods (+13 more)

### Community 7 - "Mallorca SDM"
Cohesion: 0.14
Nodes (13): 1. `data_preparation.R` (R), 2. `extract_raster_data.py` (Python), 3. `preprocessing.py` (Python), 4. `run_model.py` (Python), 5. `sampling_effort.py` (Python, opcional), Análisis y utilidades, Conjunto de especies para la evaluación, Entorno Python (+5 more)

### Community 8 - "Notas de diseño experimental"
Cohesion: 0.25
Nodes (7): Extract raster — PA no debe ser un solo punto por celda, Hold-out PAC, Notas de diseño experimental, optimal_model / transferencia PA vs PO, sampling_effort — eje X (grids vs observations), Trabajo futuro — site occupancy models (otro paper), Validación espacial PAC — clusters y leave-one-out

### Community 9 - "Current State"
Cohesion: 0.25
Nodes (7): Code, Current State, Data, Infrastructure, Manuscript, Open checklist, Scope

### Community 10 - "Market stress / idoneidad — mallorca-sdm"
Cohesion: 0.25
Nodes (7): FEASIBILITY, Market stress / idoneidad — mallorca-sdm, MARKET / TRANSFER, SCIENCE, Scope lock, Stress test actions, VERDICT

### Community 11 - "mallorca-sdm"
Cohesion: 0.29
Nodes (6): Agents, Graphify, Layout, mallorca-sdm, Paths, Scope

### Community 12 - "Vision"
Cohesion: 0.33
Nodes (5): Para quién, Qué NO es el producto (aún), Qué problema resuelve, Qué quieres conseguir, Vision

### Community 13 - "Roadmap"
Cohesion: 0.40
Nodes (4): Later, Next, Now, Roadmap

### Community 14 - "Decisions"
Cohesion: 0.50
Nodes (3): Decisions, Locked, Open (USER)

## Knowledge Gaps
- **59 isolated node(s):** `Scope`, `Agents`, `Layout`, `Paths`, `Graphify` (+54 more)
  These have ≤1 connection - possible missing edges or undocumented components. (Counts symbols only; 118 node(s) total have ≤1 connection when file, concept and rationale nodes are included.)

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `run_model()` connect `run_model.py` to `preprocessing.py`, `optimal_model.py`, `sampling_effort.py`, `models.py`?**
  _High betweenness centrality (0.034) - this node is a cross-community bridge._
- **Why does `MLP` connect `models.py` to `run_model.py`, `optimal_model.py`?**
  _High betweenness centrality (0.031) - this node is a cross-community bridge._
- **Why does `apply()` connect `optimal_model.py` to `sampling_effort.py`?**
  _High betweenness centrality (0.020) - this node is a cross-community bridge._
- **What connects `Scope`, `Agents`, `Layout` to the rest of the system?**
  _59 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `optimal_model.py` be split into smaller, more focused modules?**
  _Cohesion score 0.13538461538461538 - nodes in this community are weakly interconnected._
- **Should `models.py` be split into smaller, more focused modules?**
  _Cohesion score 0.056429232192414434 - nodes in this community are weakly interconnected._
- **Should `2. Materials and Methods` be split into smaller, more focused modules?**
  _Cohesion score 0.09090909090909091 - nodes in this community are weakly interconnected._