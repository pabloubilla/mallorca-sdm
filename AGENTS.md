# mallorca-sdm

Product Lead: **Iván Cortés-Fernández**. Colaborador: **Pablo Ubilla Pavez** (INRIA Montpellier). Portfolio hub: `~/research/portfolio_research`.  
Remote GitHub: `pabloubilla/mallorca-sdm`.

## Scope

Deep SDM multi-especie — plantas terrestres Mallorca. Dataset PA+PO con cuatro clases de observación:

| Código | Tipo | Origen (aprox.) |
|--------|------|-----------------|
| PAC | PA completo (cuadrícula) | Flora (+ filtros GBIF) |
| PAU | PA incompleto | Flora incompleta |
| POV | PO validado (punto) | Biodibal |
| POU | PO no validado (punto) | GBIF |

Benchmark: entrenar PAU/POV/POU (y mixes) → evaluar AUC por especie en **PAC**. Ver `README.md`.

## Agents

- Follow `product/` (vision, decisions, roadmap, current-state, market-stress).
- Research integrity: never invent numbers; never put author names or titles in figure artwork.
- Caveman chat unless user says otherwise (portfolio skill).
- Do not commit/push unless USER asks.
- Secrets / raw heavy data: prefer `data/` (gitignored) or NAS; never commit credentials.
- Prefer `graphify query` / `path` / `explain` before broad Grep when exploring (`graphify-out/`).

## Layout

- Root scripts: `data_preparation.R`, `extract_raster_data.py`, `preprocessing.py`, `run_model.py`, `sampling_effort.py`, `optimal_model.py`, `models.py`, …
- `data/raw/`, `data/processed/` — locales / NAS (no SoT en git)
- `output/` — checkpoints + CSV/PNG (gitignored)
- `manuscript_draft.md` — borrador; migrar a `manuscript/` LaTeX Later
- `product/` — Product Lead docs
- `references/` — PDFs literature
- `notes.md`, `checklist.md` — diseño experimental abierto

## Paths

- Local: `~/research/mallorca-sdm`
- NAS: `smb://dxp4800pro-14df.local/personal_folder/research/mallorca-sdm`

## Graphify

```bash
source /home/usuario/research/portfolio_research/.venv/bin/activate
graphify . --code-only --no-viz
# after edits:
graphify update .
```
