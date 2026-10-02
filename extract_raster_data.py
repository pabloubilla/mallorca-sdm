"""
Extract raster covariates from occurrence records.

- PAC / PAU (PA): every ~100 m pixel inside each 1×1 km Flora cell is kept as its
  own row (lon/lat at pixel centre, linked to UTMCODE1X1). No pre-aggregation.
  Enables predict-per-pixel → aggregate-after-predict downstream.
- POV / POU (PO): sample at GPS point (unchanged).
"""
from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.mask import mask
from rasterio.transform import xy

DATA_PATH = Path("data/raw/full_data_saez.gpkg")
FLORA_NET_PATH = Path("data/raw/Flora_net.gpkg")
GRID_COL = "UTMCODE1X1"
CRS = "EPSG:25831"
PA_CLASSES = ("PAC", "PAU")
PO_CLASSES = ("POV", "POU")

RASTER_FOLDERS = [
    Path("/home/usuario/research/QGIS/CTL_separated"),
    Path("/home/usuario/research/QGIS/CORINE_raster_bal/dummy/"),
    # slope.tif + elev_mallorca_25831.tif live in the same folder on this machine
    Path("/home/usuario/research/QGIS/elev_mallorca/"),
]

MAX_ROWS: int | None = None
PA_PIXELS_CSV = Path("data/raw/extracted_rasters_pa_pixels.csv")
PA_OCCURRENCES_CSV = Path("data/raw/pa_occurrences.csv")
PO_POINTS_CSV = Path("data/raw/extracted_rasters_po_points.csv")
OUTPUT_CSV = Path("data/raw/extracted_data.csv")
OUTPUT_GPKG = Path("data/raw/extracted_data.gpkg")


def _sample_points_from_raster(
    raster_path: Path, points_subset: gpd.GeoDataFrame
) -> list[float]:
    """Point sample with rasterio.sample (PO records)."""
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


def _pixel_centers_in_polygon(
    reference_raster: Path, polygon, grid_code: str
) -> pd.DataFrame:
    """
    Valid pixel centres inside a 1×1 km cell using the reference raster grid.
    Returns pixel_id, lon, lat, pixel_row, pixel_col, UTMCODE1X1.
    """
    with rasterio.open(reference_raster) as src:
        geom = gpd.GeoSeries([polygon], crs=CRS).to_crs(src.crs).iloc[0]
        try:
            out, out_transform = mask(
                src, [geom], crop=True, filled=True, nodata=np.nan
            )
        except ValueError:
            return pd.DataFrame()

        data = out[0].astype(np.float64)
        nodata = src.nodata
        if nodata is not None and not np.isnan(nodata):
            data[data == nodata] = np.nan

        valid = np.isfinite(data)
        if not valid.any():
            return pd.DataFrame()

        rows, cols = np.where(valid)
        lons, lats = zip(
            *(
                xy(out_transform, r, c, offset="center")
                for r, c in zip(rows, cols)
            )
        )
        lons = np.asarray(lons, dtype=np.float64)
        lats = np.asarray(lats, dtype=np.float64)

        # Back to project CRS for stable pixel_id / downstream join
        centres = gpd.GeoDataFrame(
            {"pixel_row": rows, "pixel_col": cols},
            geometry=gpd.points_from_xy(lons, lats, crs=src.crs),
        ).to_crs(CRS)

        out_df = pd.DataFrame(
            {
                "pixel_row": centres["pixel_row"].values,
                "pixel_col": centres["pixel_col"].values,
                GRID_COL: grid_code,
                "lon": centres.geometry.x.values,
                "lat": centres.geometry.y.values,
            }
        )
        out_df["pixel_id"] = (
            out_df[GRID_COL].astype(str)
            + "_"
            + out_df["lon"].round(3).astype(str)
            + "_"
            + out_df["lat"].round(3).astype(str)
        )
        return out_df


def _sample_rasters_at_pixels(
    pixels: pd.DataFrame,
    raster_files: list[Path],
    col_names: dict[Path, str],
) -> pd.DataFrame:
    """Add one column per raster sampled at pixel centres."""
    geom_gdf = gpd.GeoDataFrame(
        pixels,
        geometry=gpd.points_from_xy(pixels["lon"], pixels["lat"]),
        crs=CRS,
    )
    out = pixels.copy()
    for raster_path in raster_files:
        col_name = col_names[raster_path]
        out[col_name] = _sample_points_from_raster(raster_path, geom_gdf)
    return out


def collect_raster_files(folders: list[Path]) -> list[Path]:
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


def _maybe_subsample(gdf: gpd.GeoDataFrame, max_rows: int | None) -> gpd.GeoDataFrame:
    if max_rows is None:
        return gdf
    out = (
        gdf.groupby("class", group_keys=False)
        .apply(lambda x: x.sample(n=min(len(x), max_rows), random_state=42))
        .reset_index(drop=True)
    )
    print(f"Usando subconjunto de prueba: {len(out):,} filas")
    return out


def _assign_pa_grids(pa_gdf: gpd.GeoDataFrame, flora_net: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    net = flora_net[[GRID_COL, "geometry"]].copy()
    joined = gpd.sjoin(pa_gdf, net, how="left", predicate="within")
    if joined.index.duplicated().any():
        joined = joined[~joined.index.duplicated(keep="first")]
    missing = joined[GRID_COL].isna().sum()
    if missing:
        print(f"  PA warning: {missing} records outside Flora_net — dropped")
        joined = joined[joined[GRID_COL].notna()].copy()
    joined[GRID_COL] = joined[GRID_COL].astype(str)
    return joined


def _cell_class_map(pa_gdf: gpd.GeoDataFrame) -> dict[str, set[str]]:
    """Which PA classes (PAC/PAU) appear in each grid cell."""
    mapping: dict[str, set[str]] = {}
    for grid_code, grp in pa_gdf.groupby(GRID_COL):
        mapping[str(grid_code)] = set(grp["class"].unique())
    return mapping


def extract_pa_pixels(
    cells: gpd.GeoDataFrame,
    cell_classes: dict[str, set[str]],
    raster_files: list[Path],
    col_names: dict[Path, str],
) -> pd.DataFrame:
    """
    One row per raster pixel inside each active 1×1 km cell.
    Duplicated per class when a cell has both PAC and PAU records.
    """
    reference = raster_files[0]
    with rasterio.open(reference) as src:
        print(f"PA pixel grid reference: {reference.name} | res={src.res}")

    pixel_frames: list[pd.DataFrame] = []
    n_cells = len(cells)
    for i, row in enumerate(cells.itertuples(), 1):
        grid_code = str(getattr(row, GRID_COL))
        centres = _pixel_centers_in_polygon(reference, row.geometry, grid_code)
        if centres.empty:
            continue
        if i % 100 == 0 or i == n_cells:
            print(f"  cells {i}/{n_cells} | pixels so far {sum(len(f) for f in pixel_frames):,}")

        for cls in sorted(cell_classes.get(grid_code, set())):
            tagged = centres.copy()
            tagged["class"] = cls
            tagged["pixel_id"] = tagged["pixel_id"] + f"_{cls}"
            pixel_frames.append(tagged)

    if not pixel_frames:
        raise RuntimeError("No PA pixels extracted inside active cells.")

    pixels = pd.concat(pixel_frames, ignore_index=True)
    n_unique_geom = pixels.drop_duplicates(
        subset=["lon", "lat", GRID_COL, "class"]
    ).shape[0]
    print(
        f"PA pixels (with class tags): {len(pixels):,} rows | "
        f"unique (lon,lat,grid,class): {n_unique_geom:,}"
    )

    pixels = _sample_rasters_at_pixels(pixels, raster_files, col_names)
    pixels["extract_mode"] = "cell_pixel"
    return pixels


def extract_po_points(
    po_gdf: gpd.GeoDataFrame,
    raster_files: list[Path],
    col_names: dict[Path, str],
) -> tuple[pd.DataFrame, gpd.GeoDataFrame]:
    po_gdf = po_gdf.copy()
    po_gdf["_geom_key"] = po_gdf.geometry.apply(lambda g: g.wkb)

    unique_geom = (
        po_gdf.drop_duplicates(subset="_geom_key")[["_geom_key", "geometry"]]
        .reset_index(drop=True)
    )
    unique_geom = gpd.GeoDataFrame(
        unique_geom, geometry="geometry", crs=po_gdf.crs
    )
    unique_geom["geom_id"] = np.arange(len(unique_geom))
    unique_geom["lon"] = unique_geom.geometry.x
    unique_geom["lat"] = unique_geom.geometry.y

    print(
        f"PO point extract: {len(po_gdf):,} rows | "
        f"{len(unique_geom):,} unique geometries "
        f"({100 * len(unique_geom) / len(po_gdf):.1f}%)"
    )

    raster_values = unique_geom[["geom_id", "lon", "lat"]].copy()
    for raster_path in raster_files:
        col_name = col_names[raster_path]
        with rasterio.open(raster_path) as src:
            print(f"  - {col_name} | raster CRS: {src.crs}")
        n_inside = _count_inside_bounds(raster_path, unique_geom)
        print(f"      puntos dentro del extent: {n_inside}/{len(unique_geom)}")
        raster_values[col_name] = _sample_points_from_raster(
            raster_path, unique_geom
        )
        n_ok = pd.Series(raster_values[col_name]).notna().sum()
        print(f"      valores no nulos: {n_ok}/{len(raster_values)}")

    raster_values["extract_mode"] = "point"
    po_gdf = po_gdf.merge(
        unique_geom[["_geom_key", "geom_id"]], on="_geom_key", how="left"
    )
    return raster_values, po_gdf


def main() -> gpd.GeoDataFrame:
    print("Loading occurrence data")
    flora = gpd.read_file(DATA_PATH)
    if flora.crs is None:
        flora = flora.set_crs(CRS)
    print(f"CRS: {flora.crs} | species: {flora['scientificName'].nunique()}")

    flora = _maybe_subsample(flora, MAX_ROWS)

    pa = flora[flora["class"].isin(PA_CLASSES)].copy()
    po = flora[flora["class"].isin(PO_CLASSES)].copy()
    print(f"PA rows: {len(pa):,} | PO rows: {len(po):,}")

    print("Loading Flora_net (1×1 km cells)")
    flora_net = gpd.read_file(FLORA_NET_PATH)
    if flora_net.crs is None:
        flora_net = flora_net.set_crs(CRS)
    flora_net[GRID_COL] = flora_net[GRID_COL].astype(str)

    print("Raster list")
    raster_files = collect_raster_files(RASTER_FOLDERS)
    col_names = raster_column_names(raster_files)

    parts: list[pd.DataFrame] = []

    if len(pa) > 0:
        print("\n--- PA: per-pixel extract (no aggregation) ---")
        pa = _assign_pa_grids(pa, flora_net)
        cell_classes = _cell_class_map(pa)

        pa_occ = pa[
            ["class", GRID_COL, "scientificName", "eventDate"]
        ].copy()
        PA_OCCURRENCES_CSV.parent.mkdir(parents=True, exist_ok=True)
        pa_occ.to_csv(PA_OCCURRENCES_CSV, index=False)
        print(f"Saved {PA_OCCURRENCES_CSV} ({len(pa_occ):,} occurrence rows)")

        active_codes = list(cell_classes.keys())
        cells = flora_net[flora_net[GRID_COL].isin(active_codes)].copy()
        pa_pixels = extract_pa_pixels(
            cells, cell_classes, raster_files, col_names
        )
        PA_PIXELS_CSV.parent.mkdir(parents=True, exist_ok=True)
        pa_pixels.to_csv(PA_PIXELS_CSV, index=False)
        print(f"Saved {PA_PIXELS_CSV}")
        parts.append(pa_pixels)

    if len(po) > 0:
        print("\n--- PO: GPS point sample ---")
        po_rasters, po = extract_po_points(po, raster_files, col_names)
        PO_POINTS_CSV.parent.mkdir(parents=True, exist_ok=True)
        po_rasters.to_csv(PO_POINTS_CSV, index=False)
        print(f"Saved {PO_POINTS_CSV}")

        po_out = po.merge(po_rasters, on="geom_id", how="left", suffixes=("_dup", ""))
        drop_cols = [
            c
            for c in po_out.columns
            if c.endswith("_dup") or c in ("_geom_key", "geom_id")
        ]
        po_out = po_out.drop(columns=drop_cols, errors="ignore")
        po_out["extract_mode"] = "point"
        parts.append(po_out)

    if not parts:
        raise RuntimeError("No PA or PO records to extract.")

    extracted = pd.concat(parts, axis=0, ignore_index=True)

    first_col = col_names[raster_files[0]]
    n_ok = extracted[first_col].notna().sum()
    if n_ok == 0:
        raise RuntimeError(
            f"Extract vacío en {first_col}. Revisa rasters y "
            f"{PA_PIXELS_CSV} / {PO_POINTS_CSV}."
        )

    print(f"\nFilas finales: {len(extracted):,} (non-null {first_col}: {n_ok:,})")
    print(
        extracted.groupby(["class", "extract_mode"]).size()
        if "class" in extracted.columns
        else extracted["extract_mode"].value_counts()
    )

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    extracted.drop(columns="geometry", errors="ignore").to_csv(
        OUTPUT_CSV, index=False
    )
    print(f"Saved {OUTPUT_CSV}")

    out_gdf = gpd.GeoDataFrame(
        extracted,
        geometry=gpd.points_from_xy(extracted["lon"], extracted["lat"]),
        crs=CRS,
    )
    out_gdf.to_file(OUTPUT_GPKG, driver="GPKG")
    print(f"Saved {OUTPUT_GPKG}")
    print("Done.")
    return out_gdf


if __name__ == "__main__":
    main()
