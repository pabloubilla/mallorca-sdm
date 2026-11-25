# MallorcaPlant manuscript

GeoPlant-style dataset paper draft (NeurIPS Datasets & Benchmarks framing).

| File | Role |
|------|------|
| `main.tex` | **Primary draft** — MallorcaPlant |
| `main-2.tex` | Archived IDM/methods scratchpad (paper 2) |
| `references.bib` | Shared bibliography (+ `geoplant2024`) |
| `figures/` | Map, mix AUC, sampling-effort curves |
| `neurips_2024.sty` | Official NeurIPS 2024 style (optional; needs full TeX Live) |

## Compile

```bash
cd manuscript/mallorca-sdm-manuscript
# if system TEXMFVAR is not writable:
export TEXMFVAR="$(pwd)/.texmf-var" TEXMFCACHE="$(pwd)/.texmf-var"
pdflatex main
pdflatex main
```

`main.tex` currently uses a self-contained `article` preamble so it builds without Times/`environ`. On Overleaf or full TeX Live, you can switch the preamble to `\usepackage[preprint]{neurips_2024}` (sty already vendored).
