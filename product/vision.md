# Vision

## Qué problema resuelve

Los Deep SDM multi-especie mezclan atlas (PA) y citizen science (PO) sin protocolo claro: qué fuente, qué resolución espacial, qué loss, y cómo escala el rendimiento con el esfuerzo. Mallorca permite un laboratorio controlado con Flora (PA completo/incompleto), Biodibal (PO validado) y GBIF (PO no validado).

## Para quién

- Comunidad científica (SDM, deep learning, biogeografía, citizen science)
- Curadores / gestores de atlas y Biodibal (interpretación de sesgos de muestreo)
- Portfolio UIB / coautores del dataset Mallorca

## Qué quieres conseguir

1. **Paper:** benchmark transferible — train PAU/POV/POU → eval PAC; curvas de esfuerzo; mixes de fuentes; point vs grouped.
2. **Dataset:** PA+PO grande Mallorca, taxonomía unificada (GBIF backbone), clases PAC/PAU/POV/POU explícitas.
3. **Método honesto:** loss por tipo de dato (BCE vs DeepMaxEnt); CV espacial en PAC; no tratar ceros PO como ausencias.

## Qué NO es el producto (aún)

- No es SaaS ni app de campo.
- No es reevaluación IUCN ni mapa de gestión operativa por especie.
- No es el paper de occupancy/detection (`notes.md` → Later).
- No sustituye muestreo Flora ni validación Biodibal.
