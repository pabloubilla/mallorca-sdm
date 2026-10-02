# Notas de diseño experimental

## optimal_model / transferencia PA vs PO

- Los 0 en PO no son ausencias; mezclar PAU+PO con BCE global trata ceros PO como ausencias reales → incorrecto.
- Comparar brazos: (A) todo PO con MaxEnt; (C) loss por fila según fuente (PAU→BCE, POV/POU→solo presencias); (D) PAC como oráculo.
- Mezclas 0/50/100 % útiles para cantidad; primero fijar **estrategia de loss**, luego esfuerzo.
- Entrenamiento comparado en **celdas** (`_grouped` para PO); evaluación siempre en PAC (PA completo).

## Hold-out PAC

- Partir PAC en `PAC_train` / `PAC_test` (idealmente bloques espaciales, no aleatorio i.i.d.).
- Métricas **solo en PAC_test**. Excluir esas celdas (y puntos PO dentro) del train en PAU/POV/POU.
- Brazos: T0 sin PAC; T1 solo PAC_train; T2 PAC_train + PO informado; comparar si entrenar con PAC mejora vs transferencia.
- ~366 celdas PAC → considerar k-fold espacial además de un split fijo.

## Validación espacial PAC — clusters y leave-one-out

- AUC en todo PAC sin CV espacial **sobreestima** el rendimiento: celdas vecinas comparten bioclima y especies (autocorrelación espacial).
- **Clusterizar PAC por coordenadas** (`lon`/`lat` o centroide de `UTMCODE1X1`): p. ej. k-means, hierarchical, o bloques espaciales regulares en UTM.
- En cada iteración: dejar **un cluster fuera** (leave-one-cluster-out, LOCO) o un bloque fuera (spatial LOO / k-fold por cluster); entrenar en el resto de PAC (y/o otras fuentes) y evaluar AUC **solo en las celdas del cluster hold-out**.
- Repetir para todos los clusters → media ± IC de AUC por especie y global; estimación más honesta que un único split aleatorio por celda.
- Al entrenar con PAU/POV/POU: excluir del train las celdas (y puntos PO dentro) que caigan en el cluster PAC de test de esa iteración.
- Con ~366 celdas y pocos clusters, el número de folds y el tamaño de cada hold-out hay que equilibrarlos (p. ej. 5–10 clusters, no LOO celda a celda salvo análisis exploratorio).

## sampling_effort — eje X (grids vs observations)

- Grids, filas de matriz (`n_train_rows`) y observaciones crudas (`n_obs`) no son intercambiables sin redefinir el subsampling.
- PAU/PAC: unidad natural = **celda 1×1 km**; PO: **registros / sitios GPS**.
- No eliminar matrices grid para PA; para PO con eje en observations usar matriz punto + MaxEnt.
- Paneles separados recomendados: PAU por celdas; POV/POU por observations o sitios; grouped PO para comparación espacial con PA.
- Si un solo eje `n_obs` para todas las fuentes: PAU debe agregarse por celda tras muestrear registros y documentar que el presupuesto es “registros de campo”, no mismo diseño espacial.

## Extract raster — PA no debe ser un solo punto por celda

- Rasters a **~100 m**; Flora/PA en cuadrícula **1×1 km** (`UTMCODE1X1`).
- Extraer covariables en un único punto por registro PA (centroide o coordenada del registro) **no representa la celda**: mezcla resoluciones y submuestrea el raster dentro de la celda.
- Para PA (PAC/PAU): **un registro por píxel raster (~100 m)** dentro de cada celda 1×1 km (`pixel_id`, `lon`/`lat` del centro del píxel, `UTMCODE1X1`). Sin media en el extract.
- Etiquetas de especie a nivel de **celda** (broadcast a todos los píxeles de esa celda en `preprocessing.py`).
- Inferencia: predecir en cada píxel → **agregar a celda** (media de probabilidades) → AUC en PAC (`run_model._aggregate_pixel_scores_to_cells`).
- PO sigue en GPS punto. Re-ejecutar `extract_raster_data.py` + `preprocessing.py` tras cambiar rasters.

## Trabajo futuro — site occupancy models (otro paper)

- Experimentar con **modelos de ocupación por sitio** (occupancy + detection) para estimar **probabilidad de detección** por grupo taxonómico y por **método/fuente** (Flora atlas, Biodibal, GBIF).
- Encaja con la distinción PA vs PO y con PAU incompleto vs PAC completo: separar ocupación verdadera de probabilidad de registrar presencia/ausencia según el protocolo.
- Requiere diseño explícito de grupos taxonómicos y estratificación por método de muestreo; no sustituye el benchmark deep SDM actual pero complementa la interpretación de sesgos de detección.
