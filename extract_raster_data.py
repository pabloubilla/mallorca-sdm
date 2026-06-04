#### Extract raster data from Flora data
from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio

DATA_PATH = "./data/raw/full_data_clean_saez.gpkg"
# Una o varias carpetas con archivos .tif (recursivo solo en el nivel de cada carpeta)
RASTER_FOLDERS = [
    Path("/Volumes/Crucial_SSD/QGIS/Climate_UIB/CTL_separated"),
    Path("/Volumes/Crucial_SSD/QGIS/CORINE_EU/CORINE_raster_bal/dummy"),
    Path("/Volumes/Crucial_SSD/QGIS/World_ELE_GISdata_GlobalSolarAtlas-v2_GEOTIFF/MDT05_Mallorca/Pendiente"),
    Path("/Volumes/Crucial_SSD/QGIS/World_ELE_GISdata_GlobalSolarAtlas-v2_GEOTIFF/elev_mallorca/"),
]

# Prueba con subconjunto (poner None para procesar todo)
MAX_ROWS = None
UNIQUE_EXTRACT_CSV = "./data/raw/extracted_rasters_unique_geom.csv"


def _sample_points_from_raster(
    raster_path: Path, points_subset: gpd.GeoDataFrame
) -> list[float]:
    """Muestreo con rasterio.sample (equivalente fiable a extract de R)."""
    with rasterio.open(raster_path) as src:
        pts = points_subset.to_crs(src.crs)
        coords = [(point.x, point.y) for point in pts.geometry]

        values: list[float] = []
        nodata = src.nodata
        for val in src.sample(coords):
            v = float(val[0])
            if nodata is not None and not np.isnan(nodata) and v == nodata:
                v = np.nan
            elif np.isnan(v):
                v = np.nan
            values.append(v)

        return values


def collect_raster_files(folders: list[Path]) -> list[Path]:
    """Lista todos los .tif de cada carpeta (ordenados por ruta)."""
    raster_files: list[Path] = []
    for folder in folders:
        folder = Path(folder)
        if not folder.is_dir():
            raise FileNotFoundError(f"Carpeta no encontrada: {folder}")
        found = sorted(folder.glob("*.tif"))
        if not found:
            print(f"  (aviso) sin .tif en {folder}")
        else:
            print(f"  {len(found)} .tif en {folder}")
        raster_files.extend(found)
    if not raster_files:
        folders_str = ", ".join(str(f) for f in folders)
        raise FileNotFoundError(f"No hay .tif en ninguna carpeta: {folders_str}")
    return raster_files


def raster_column_names(raster_files: list[Path]) -> dict[Path, str]:
    """Nombre de columna por raster; prefijo de carpeta si el stem se repite."""
    stems = [p.stem for p in raster_files]
    stem_counts = pd.Series(stems).value_counts()
    names: dict[Path, str] = {}
    for path in raster_files:
        if stem_counts[path.stem] > 1:
            names[path] = f"{path.parent.name}_{path.stem}"
        else:
            names[path] = path.stem
    return names


def _count_inside_bounds(raster_path: Path, points_subset: gpd.GeoDataFrame) -> int:
    with rasterio.open(raster_path) as src:
        pts = points_subset.to_crs(src.crs)
        b = src.bounds
        inside = [
            b.left <= p.x <= b.right and b.bottom <= p.y <= b.top
            for p in pts.geometry
        ]
    return int(sum(inside))


print("Loading Flora data")
Flora_data = gpd.read_file(DATA_PATH)
print("Number of unique species is", len(Flora_data['scientificName'].unique()))
print(f"CRS puntos: {Flora_data.crs}")

if MAX_ROWS is not None:
    Flora_data = (

    Flora_data

    .groupby("class", group_keys=False)

    .apply(lambda x: x.sample(n=min(len(x), MAX_ROWS), random_state=42))

    .reset_index(drop=True)

)
    print(f"Usando subconjunto de prueba: {len(Flora_data):,} filas")

# Geometrías únicas (muchas filas comparten punto)
Flora_data["_geom_key"] = Flora_data.geometry.apply(lambda g: g.wkb)

unique_geom = (
    Flora_data.drop_duplicates(subset="_geom_key")[["_geom_key", "geometry"]]
    .reset_index(drop=True)
)
unique_geom = gpd.GeoDataFrame(unique_geom, geometry="geometry", crs=Flora_data.crs)
unique_geom["geom_id"] = np.arange(len(unique_geom))
# Coordenadas en EPSG:25831 (mismo CRS que puntos y rasters)
unique_geom["lon"] = unique_geom.geometry.x
unique_geom["lat"] = unique_geom.geometry.y

print(
    f"Filas totales: {len(Flora_data):,} | "
    f"Geometrías únicas: {len(unique_geom):,} "
    f"({100 * len(unique_geom) / len(Flora_data):.1f}%)"
)

print("Getting list of raster files")
raster_files = collect_raster_files(RASTER_FOLDERS)
col_names = raster_column_names(raster_files)

raster_values = unique_geom[["geom_id", "lon", "lat"]].copy()

print(f"Extracting with rasterio.sample from {len(raster_files)} rasters")
for raster_path in raster_files:
    col_name = col_names[raster_path]
    with rasterio.open(raster_path) as src:
        print(f"  - {col_name} | raster CRS: {src.crs}")
    n_inside = _count_inside_bounds(raster_path, unique_geom)
    print(f"      puntos dentro del extent: {n_inside}/{len(unique_geom)}")

    raster_values[col_name] = _sample_points_from_raster(raster_path, unique_geom)
    n_ok = pd.Series(raster_values[col_name]).notna().sum()
    print(f"      valores no nulos: {n_ok}/{len(raster_values)}")

# Guardar SIEMPRE el extract antes del merge
print(f"Saving unique-geometry extract backup -> {UNIQUE_EXTRACT_CSV}")
raster_values.to_csv(UNIQUE_EXTRACT_CSV, index=False)

first_col = col_names[raster_files[0]]
if raster_values[first_col].notna().sum() == 0:
    raise RuntimeError(
        f"Extract vacío en {first_col} aunque CRS coincide. "
        "Revisa extent del raster y coordenadas en {UNIQUE_EXTRACT_CSV}."
    )

# Join de vuelta a todas las filas
Flora_with_geom_id = Flora_data.merge(
    unique_geom[["geom_id", "_geom_key"]],
    on="_geom_key",
    how="left",
)

extracted_data = Flora_with_geom_id.merge(
    raster_values,
    on="geom_id",
    how="left",
).drop(columns=["_geom_key", "geom_id"])
# lon, lat en EPSG:25831 (metros, UTM zona 31N)

print(extracted_data.head())
print(f"Filas finales: {len(extracted_data):,}")

print("Saving extracted data to csv")
extracted_data.drop(columns="geometry", errors="ignore").to_csv(
     "data/raw/extracted_data.csv", index=False
)

print("Saving extracted data to gpkg")
extracted_data.to_file("data/raw/extracted_data.gpkg", driver="GPKG")

print("Done.")
