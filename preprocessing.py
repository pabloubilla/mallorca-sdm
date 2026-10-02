"""Processed matrices: PA by raster pixel; PO by point + grid-mean grouped (_grouped)."""
import os

import geopandas as gpd
import numpy as np
import pandas as pd

from covariables import COVARS_LIST

FLORA_NET_PATH = "data/raw/Flora_net.gpkg"
EXTRACTED_CSV = "data/raw/extracted_data.csv"
PA_OCCURRENCES_CSV = "data/raw/pa_occurrences.csv"
GRID_COL = "UTMCODE1X1"
GROUPED_SUFFIX = "_grouped"
PIXEL_ID_COL = "pixel_id"

PA_CLASSES = ("PAC", "PAU")
PO_CLASSES = ("POV", "POU")

MIN_SPECIES_LIST_PAC_GRIDS = 1
MIN_TRAIN_OCC_GRID = 1
MIN_TRAIN_OCC_POINT = 2


def file_suffix(source: str, grouped: bool = False) -> str:
    if grouped and is_po_source(source):
        return GROUPED_SUFFIX
    return ""


def species_matrix_path(source: str, grouped: bool = False) -> str:
    return f"data/processed/species_matrix_{source}{file_suffix(source, grouped)}.csv"


def covariates_path(source: str, grouped: bool = False) -> str:
    return f"data/processed/covariates_{source}{file_suffix(source, grouped)}.csv"


def grid_cells_path(source: str, grouped: bool = False) -> str:
    return f"data/processed/grid_cells_{source}{file_suffix(source, grouped)}.csv"


def is_po_source(source: str) -> bool:
    return source in PO_CLASSES


def assign_cod10_grid(df: pd.DataFrame) -> pd.DataFrame:
    gdf = gpd.GeoDataFrame(
        df,
        geometry=gpd.points_from_xy(df["lon"], df["lat"]),
        crs="EPSG:25831",
    )
    # Extract concat leaves empty UTMCODE1X1 on PO; drop so sjoin keeps GRID_COL name.
    if GRID_COL in gdf.columns:
        gdf = gdf.drop(columns=[GRID_COL])
    net = gpd.read_file(FLORA_NET_PATH)[[GRID_COL, "geometry"]]
    joined = gpd.sjoin(gdf, net, how="left", predicate="within")
    if joined.index.duplicated().any():
        joined = joined[~joined.index.duplicated(keep="first")]
    missing = joined[GRID_COL].isna().sum()
    if missing:
        print(f"  warning: {missing} rows outside Flora_net — dropped")
        joined = joined[joined[GRID_COL].notna()].copy()
    return pd.DataFrame(joined.drop(columns="geometry"))


def _point_site_id(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["site_id"] = out["lon"].astype(str) + "_" + out["lat"].astype(str)
    return out


def _species_matrix(
    subset: pd.DataFrame, id_col: str, species_list: list[str]
) -> pd.DataFrame:
    return (
        subset.groupby([id_col, "scientificName"])
        .size()
        .unstack(fill_value=0)
        .clip(upper=1)
        .astype(np.int8)
        .reindex(columns=species_list, fill_value=0)
    )


def _cell_species_matrix(
    occurrences: pd.DataFrame, species_list: list[str]
) -> pd.DataFrame:
    return _species_matrix(occurrences, GRID_COL, species_list)


def _broadcast_cell_labels_to_pixels(
    pixels: pd.DataFrame, cell_species: pd.DataFrame
) -> pd.DataFrame:
    """Attach cell-level 0/1 species vector to each raster pixel in that cell."""
    cell_species = cell_species.copy()
    cell_species.index = cell_species.index.astype(str)
    pixels = pixels.copy()
    pixels[GRID_COL] = pixels[GRID_COL].astype(str)

    labels = cell_species.reindex(pixels[GRID_COL].values)
    labels.index = pixels[PIXEL_ID_COL].values
    labels.index.name = "site_id"
    return labels


def _save_pa_pixel_products(
    pixels: pd.DataFrame,
    occurrences: pd.DataFrame,
    cls: str,
    species_list: list[str],
    covar_cols: list[str],
) -> tuple[int, int]:
    """
    One row per ~100 m pixel; species labels from 1×1 km cell (same for all
    pixels in the cell). Predict per pixel, aggregate to cell after predict.
    """
    pixels = pixels[pixels["class"] == cls].copy()
    if pixels.empty:
        print(f"No pixel rows for class {cls}")
        return 0, 0

    occ = occurrences[occurrences["class"] == cls]
    cell_species = _cell_species_matrix(occ, species_list)
    species_matrix = _broadcast_cell_labels_to_pixels(pixels, cell_species)
    species_matrix.to_csv(species_matrix_path(cls, grouped=False))

    cov = pixels.set_index(PIXEL_ID_COL)[covar_cols].copy()
    cov.index.name = "site_id"
    cov[GRID_COL] = pixels.set_index(PIXEL_ID_COL)[GRID_COL]
    cov.to_csv(covariates_path(cls, grouped=False))

    pd.DataFrame(
        {GRID_COL: pixels[GRID_COL].astype(str).unique()}
    ).to_csv(grid_cells_path(cls, grouped=False), index=False)

    n_pixels = len(species_matrix)
    n_cells = pixels[GRID_COL].nunique()
    print(
        f"\n[{cls}]  PIXEL  {n_pixels:,} pixels ({n_cells} cells) × "
        f"{species_matrix.shape[1]} species  |  {len(occ):,} occurrence rows"
    )
    return n_cells, n_pixels


def _save_grid_products(
    subset: pd.DataFrame,
    cls: str,
    species_list: list[str],
    covar_cols: list[str],
    suffix: str = "",
) -> tuple[int, int]:
    """One row per grid cell (PO grouped mode)."""
    grid_id = GRID_COL
    species_matrix = _species_matrix(subset, grid_id, species_list)
    species_matrix.index.name = "site_id"
    species_matrix.to_csv(species_matrix_path(cls, grouped=bool(suffix)))

    cov = subset.groupby(grid_id)[covar_cols].mean()
    cov.index.name = "site_id"
    cov.to_csv(covariates_path(cls, grouped=bool(suffix)))

    species_matrix.index.to_frame(index=False).to_csv(
        grid_cells_path(cls, grouped=bool(suffix)), index=False
    )

    n_grids = len(species_matrix)
    label = "GROUPED" if suffix else "GRID"
    print(
        f"\n[{cls}]  {label}  {n_grids} cells × {species_matrix.shape[1]} species"
        f"  |  {len(subset):,} obs"
    )
    return n_grids, len(subset)


def _save_point_products(
    subset: pd.DataFrame,
    cls: str,
    species_list: list[str],
    covar_cols: list[str],
) -> None:
    subset = _point_site_id(subset)
    species_matrix = _species_matrix(subset, "site_id", species_list)
    species_matrix.to_csv(species_matrix_path(cls, grouped=False))

    cov = subset.groupby("site_id")[covar_cols].first()
    cov[GRID_COL] = subset.groupby("site_id")[GRID_COL].first()
    cov.to_csv(covariates_path(cls, grouped=False))

    pd.DataFrame({GRID_COL: subset[GRID_COL].astype(str).unique()}).to_csv(
        grid_cells_path(cls, grouped=False), index=False
    )

    print(
        f"\n[{cls}]  POINT {len(species_matrix):,} sites × "
        f"{species_matrix.shape[1]} species  |  {len(subset):,} obs"
    )
    _save_grid_products(
        subset, cls, species_list, covar_cols, suffix=GROUPED_SUFFIX
    )


def main():
    df_raw = pd.read_csv(EXTRACTED_CSV)
    covar_cols = [c for c in COVARS_LIST if c in df_raw.columns]
    missing_covars = [c for c in ("lon", "lat") if c not in covar_cols]
    if missing_covars:
        raise ValueError(f"extracted_data.csv missing spatial columns: {missing_covars}")

    # Drop rows with any covariate NA (slope + CORINE gaps) so StandardScaler stays finite.
    df = df_raw.dropna(subset=covar_cols).copy()
    print(f"rows dropped (covariate NA): {df_raw.shape[0] - df.shape[0]}")

    pa_pixels = df[df.get("extract_mode", "") == "cell_pixel"].copy()
    po_rows = df[df.get("extract_mode", "") != "cell_pixel"].copy()

    if not po_rows.empty:
        print(f"Assigning {GRID_COL} to PO rows from Flora_net …")
        po_rows = assign_cod10_grid(po_rows)

    if not os.path.isfile(PA_OCCURRENCES_CSV):
        raise FileNotFoundError(
            f"Missing {PA_OCCURRENCES_CSV}. Re-run extract_raster_data.py."
        )
    pa_occ = pd.read_csv(PA_OCCURRENCES_CSV)
    pa_occ[GRID_COL] = pa_occ[GRID_COL].astype(str)

    os.makedirs("data/processed", exist_ok=True)

    pac_occ = pa_occ[pa_occ["class"] == "PAC"]
    grids_per_sp = pac_occ.groupby("scientificName")[GRID_COL].nunique()
    species_list = grids_per_sp[
        grids_per_sp >= MIN_SPECIES_LIST_PAC_GRIDS
    ].index.tolist()
    print(
        f"species catalog: {len(species_list)} "
        f"(PAC ≥{MIN_SPECIES_LIST_PAC_GRIDS} cells)"
    )

    for cls in PA_CLASSES:
        if pa_pixels.empty:
            print(f"No PA pixel data for {cls}")
            continue
        _save_pa_pixel_products(
            pa_pixels, pa_occ, cls, species_list, covar_cols
        )

    for cls in PO_CLASSES:
        subset = po_rows[po_rows["class"] == cls].copy()
        if subset.empty:
            print(f"No data for class {cls}")
            continue
        _save_point_products(subset, cls, species_list, covar_cols)


if __name__ == "__main__":
    main()
