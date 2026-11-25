# Vision

## Qué problema resuelve

Los Deep SDM multi-especie mezclan atlas (PA) y citizen science (PO) sin protocolo claro: qué fuente, qué resolución espacial, qué loss, y cómo escala el rendimiento con el esfuerzo. La dicotomía binaria PO vs PA oculta **niveles** de calidad de observación. Mallorca permite un laboratorio controlado con Flora (PA completo/incompleto), Biodibal (PO validado) y GBIF (PO no validado).

## Para quién

- Comunidad científica (SDM, deep learning, biogeografía, citizen science, IDM)
- Curadores / gestores de atlas y Biodibal (interpretación de sesgos de muestreo)
- Portfolio UIB / coautores

## Qué quieres conseguir

1. **Paper canónico — Beyond PO and PA:** evidenciar que hay **niveles** de observación (PAC / PAU / POV / POU), no solo presencia-ausencia vs presencia-solo. RQ: (1) rareza / sesgo por fuente; (2) valor de invertir en PA incompleto; (3) mixes óptimos para inventario **y** para transferencia SDM a PAC.
2. **Base curada (insumo, no el paper):** PA+PO Mallorca, taxonomía GBIF, filtro flora Mallorca ∪ listado oficial alóctonos; resto excluido como ruido.
3. **Método honesto:** loss por tipo de dato (BCE vs DeepMaxEnt); comparación justa a grano de celda (PO grouped); CV espacial en PAC (Next); no tratar ceros PO como ausencias.

## Qué NO es el producto (aún)

- **No es un paper-dataset** (GeoPlant / MallorcaPlant descartado como framing principal).
- No es SaaS ni app de campo.
- No es reevaluación IUCN ni mapa de gestión operativa por especie.
- No es el paper de occupancy/detection / second spatial field (`notes.md` → Later).
- No sustituye muestreo Flora ni validación Biodibal.
