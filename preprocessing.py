"""Processed matrices: PA by grid; PO by point (default) + grid-mean grouped (_grouped)."""
import os

import geopandas as gpd
import numpy as np
import pandas as pd

from covariables import COVARS_LIST

FLORA_NET_PATH = "data/raw/Flora_net.gpkg"
EXTRACTED_CSV = "data/raw/extracted_data.csv"
GRID_COL = "UTMCODE1X1"
GROUPED_SUFFIX = "_grouped"

PA_CLASSES = ("PAC", "PAU")
PO_CLASSES = ("POV", "POU")

MIN_SPECIES_LIST_PAC_GRIDS = 1 # minimum number of cells in PAC to include a species in the catalog
MIN_TRAIN_OCC_GRID = 1 # minimum number of cells in the training set to include a species in the training set
MIN_TRAIN_OCC_POINT = 2 # minimum number of points in the training set to include a species in the training set


def file_suffix(source: str, grouped: bool = False) -> str:
    """PO only: '_grouped' files. PA (PAC/PAU) always use the grid files without suffix."""
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


def _save_grid_products(
    subset: pd.DataFrame,
    cls: str,
    species_list: list[str],
    bioclim_cols: list[str],
    suffix: str = "",
) -> tuple[int, int]:
    """One row per grid cell; suffix '' or '_grouped'."""
    grid_id = GRID_COL
    species_matrix = _species_matrix(subset, grid_id, species_list)
    species_matrix.index.name = "site_id"
    species_matrix.to_csv(species_matrix_path(cls, grouped=bool(suffix)))

    cov = subset.groupby(grid_id)[bioclim_cols].mean()
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
    bioclim_cols: list[str],
) -> None:
    """Point-level files (no suffix) + grouped grid aggregate for PO."""
    subset = _point_site_id(subset)
    species_matrix = _species_matrix(subset, "site_id", species_list)
    species_matrix.to_csv(species_matrix_path(cls, grouped=False))

    cov = subset.groupby("site_id")[bioclim_cols].first()
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
        subset, cls, species_list, bioclim_cols, suffix=GROUPED_SUFFIX
    )


def main():
    df_raw = pd.read_csv(EXTRACTED_CSV)
    df = df_raw[df_raw["slope"].notna()].copy()
    print("rows dropped (slope NA):", df_raw.shape[0] - df.shape[0])

    print(f"Assigning {GRID_COL} from Flora_net …")
    df = assign_cod10_grid(df)

    bioclim_cols = [c for c in COVARS_LIST if c in df.columns]
    os.makedirs("data/processed", exist_ok=True)

    pac = df[df["class"] == "PAC"]
    grids_per_sp = pac.groupby("scientificName")[GRID_COL].nunique()
    species_list = grids_per_sp[
        grids_per_sp >= MIN_SPECIES_LIST_PAC_GRIDS
    ].index.tolist()
    print(
        f"species catalog: {len(species_list)} "
        f"(PAC ≥{MIN_SPECIES_LIST_PAC_GRIDS} cells)"
    )

    for cls in PA_CLASSES + PO_CLASSES:
        subset = df[df["class"] == cls].copy()
        if subset.empty:
            print(f"No data for class {cls}")
            continue
        if cls in PA_CLASSES:
            _save_grid_products(subset, cls, species_list, bioclim_cols)
        else:
            _save_point_products(subset, cls, species_list, bioclim_cols)


if __name__ == "__main__":
    main()
