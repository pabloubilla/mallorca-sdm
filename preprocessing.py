"""Processed matrices: PA by COD10X10 grid; PO by point (GPS), effort still in grids."""
import os

import geopandas as gpd
import numpy as np
import pandas as pd

from covariables import COVARS_LIST

FLORA_NET_PATH = "data/raw/Flora_net.gpkg"
EXTRACTED_CSV = "data/raw/extracted_data.csv"
GRID_COL = "UTMCODE1X1"

PA_CLASSES = ("PAC", "PAU")
PO_CLASSES = ("POV", "POU")

# Species catalog: min presences on PAC *grid cells*
MIN_SPECIES_LIST_PAC_GRIDS = 3
# Training filters (applied in run_model on the loaded subset)
MIN_TRAIN_OCC_GRID = 10
MIN_TRAIN_OCC_POINT = 10


def is_po_source(source: str) -> bool:
    return source in PO_CLASSES


def is_pa_source(source: str) -> bool:
    return source in PA_CLASSES or source == "PAC"


def grid_cells_path(source: str) -> str:
    return f"data/processed/grid_cells_{source}.csv"


def assign_cod10_grid(df: pd.DataFrame) -> pd.DataFrame:
    """Spatial join lon/lat (EPSG:25831) to Flora_net; one COD10X10 per row."""
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
) -> tuple[int, int]:
    """One row per COD10X10; returns (n_grids, n_obs)."""
    grid_id = GRID_COL
    species_matrix = _species_matrix(subset, grid_id, species_list)
    species_matrix.index.name = "site_id"
    species_matrix.to_csv(f"data/processed/species_matrix_{cls}.csv")

    cov = subset.groupby(grid_id)[bioclim_cols].mean()
    cov.index.name = "site_id"
    cov.to_csv(f"data/processed/covariates_{cls}.csv")

    grids = species_matrix.index.to_frame(index=False)
    grids.to_csv(grid_cells_path(cls), index=False)

    n_obs = len(subset)
    n_grids = len(species_matrix)
    print(
        f"\n[{cls}]  GRID  {n_grids} cells × {species_matrix.shape[1]} species"
        f"  |  {n_obs:,} obs  ({n_obs / max(n_grids, 1):.1f} obs/cell)"
    )
    return n_grids, n_obs


def _save_point_products(
    subset: pd.DataFrame,
    cls: str,
    species_list: list[str],
    bioclim_cols: list[str],
) -> tuple[int, int]:
    """One row per GPS point; grid index for sampling_effort."""
    subset = _point_site_id(subset)
    n_obs = len(subset)

    species_matrix = _species_matrix(subset, "site_id", species_list)
    species_matrix.to_csv(f"data/processed/species_matrix_{cls}.csv")

    cov = subset.groupby("site_id")[bioclim_cols].first()
    cov[GRID_COL] = subset.groupby("site_id")[GRID_COL].first()
    cov.to_csv(f"data/processed/covariates_{cls}.csv")

    grid_sub = subset.copy()
    grid_sub["site_id"] = grid_sub[GRID_COL].astype(str)
    grid_matrix = _species_matrix(grid_sub, "site_id", species_list)
    grid_matrix.to_csv(f"data/processed/species_matrix_{cls}_grid.csv")

    grids = pd.DataFrame({GRID_COL: subset[GRID_COL].astype(str).unique()})
    grids.to_csv(grid_cells_path(cls), index=False)

    n_grids = len(grids)
    n_points = len(species_matrix)
    print(
        f"\n[{cls}]  POINT {n_points:,} sites × {species_matrix.shape[1]} species"
        f"  |  {n_obs:,} obs  |  {n_grids} grids for sampling effort"
    )
    return n_grids, n_obs


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
        f"(PAC ≥{MIN_SPECIES_LIST_PAC_GRIDS} grid cells)"
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
