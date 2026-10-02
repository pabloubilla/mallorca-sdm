# Market stress / idoneidad — mallorca-sdm

Date: 2026-09-25  
Owner: Product Lead (TBD) · Portfolio Lead viability frame

## Scope lock

**Mallorca · plantas terrestres · Deep SDM multi-especie.**  
Fuentes: Flora (PAC/PAU), Biodibal (POV), GBIF (POU).  
Eval: transferencia a PAC. No SaaS.

## SCIENCE

- Niche: efecto del **tipo y calidad de dato** (PA completo/incompleto × PO validado/no) + esfuerzo + resolución (point vs grid) en Deep SDM.
- Overlap: MaxEnt / deep SDM / citizen-science bias literature — diferenciador = cuatro clases explícitas + eval en atlas completo controlado + loss por fuente.
- Publishable si CV espacial honesta + cobertura isla + métricas reproducibles (no overclaim del draft).

## MARKET / TRANSFER

- Users: investigadores SDM; curadores atlas/Biodibal.
- Who pays: ciencia / proyectos UIB — no claim de producto comercial.
- Market stress: **low** (paper); **high** si se vende como herramienta operativa sin validación espacial.

## FEASIBILITY

- Código + datos locales ya existen; pipeline documentado.
- Blockers: cobertura extract, CV espacial, Product Lead/coautores, LaTeX.

## VERDICT

**go** (paper metodológico) · **hold** (transferencia operativa / segundo paper occupancy)

## Stress test actions

1. Re-run effort + optimal tras ampliar celdas; citar solo `output/`.
2. Implementar LOCO PAC.
3. `literature-gap` + journal scout cuando Product Lead cerrado.
