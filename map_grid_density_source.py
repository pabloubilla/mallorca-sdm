"""2×2 map: observation density per UTMCODE1X1 cell for PAC, PAU, POV, POU.

Uses ``full_data_saez.gpkg`` plus Biel FLORA GBO Mallorca (PAC/PAU from
cuadrículas GBO), so the figure reflects Biel even before re-running extract.
"""
from __future__ import annotations

import re
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from plot_style import apply as apply_plot_style
from preprocessing import FLORA_NET_PATH, GRID_COL

SOURCES = ["PAC", "PAU", "POV", "POU"]
PANEL_LETTERS = ["A", "B", "C", "D"]
FULL_DATA_GPKG = Path("data/raw/full_data_saez.gpkg")
BIEL_QUAD = Path("data/raw/Biel/20260620_Quadrícules_GBO.xlsx")
BIEL_FLORA = Path("data/raw/Biel/20260823_FLORA_GBO_Reduida.xlsx")
OUTPUT_PNG = Path("output/grid_density_by_source.png")
OUTPUT_PDF = Path("output/grid_density_by_source.pdf")

# Flora_net spans Pitiusas–Menorca. Mallorca mainland ends ~541.5 km E;
# next cells jump to ~567.5 km (channel / Menorca). Clip before plot.
MALLORCA_X_MIN = 430_000.0
MALLORCA_X_MAX = 545_000.0

_FLORA_UTM_RE = re.compile(
    r"^31S\s+([A-Z]{2})\s+(\d+)\s+(\d+)(?:\s+([A-Da-d]))?$",
    re.IGNORECASE,
)


def _clip_mallorca(net: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Keep only 1 km cells on Mallorca (drop Pitiusas + Menorca)."""
    cent_x = net.geometry.centroid.x
    out = net.loc[(cent_x >= MALLORCA_X_MIN) & (cent_x <= MALLORCA_X_MAX)].copy()
    print(f"Flora_net Mallorca clip: {len(out):,} / {len(net):,} cells")
    return out


def utm_quad_to_code(utm: str) -> str | None:
    """'DD 459 4373' → 'DD5973' (same rule as data_preparation.R)."""
    if utm is None or (isinstance(utm, float) and np.isnan(utm)):
        return None
    parts = str(utm).strip().split()
    if len(parts) < 3:
        return None
    return f"{parts[0]}{parts[1][-2:]}{parts[2][-2:]}"


def parse_flora_utm_code(utm_1x1: str) -> str | None:
    """'31S DE 486 4404 A' → 'DE8644'."""
    if utm_1x1 is None or (isinstance(utm_1x1, float) and np.isnan(utm_1x1)):
        return None
    m = _FLORA_UTM_RE.match(str(utm_1x1).strip())
    if not m:
        return None
    letters, e, n, _quad = m.groups()
    return f"{letters.upper()}{e[-2:]}{n[-2:]}"


def _mallorca_limits(
    net: gpd.GeoDataFrame, active_codes: pd.Index
) -> tuple[float, float, float, float]:
    """Extent of grid cells with data, padded (~Mallorca focus)."""
    sub = net[net[GRID_COL].isin(active_codes)]
    if sub.empty:
        return tuple(net.total_bounds)
    xmin, ymin, xmax, ymax = sub.total_bounds
    pad_x = (xmax - xmin) * 0.05 or 5000
    pad_y = (ymax - ymin) * 0.05 or 5000
    return xmin - pad_x, ymin - pad_y, xmax + pad_x, ymax + pad_y


def _load_gpkg_grid_counts(net: gpd.GeoDataFrame) -> pd.DataFrame:
    gdf = gpd.read_file(FULL_DATA_GPKG)[["class", "geometry"]]
    if gdf.crs is None:
        gdf = gdf.set_crs(25831)
    elif gdf.crs.to_epsg() != 25831:
        gdf = gdf.to_crs(25831)
    joined = gpd.sjoin(gdf, net, how="inner", predicate="within")
    if joined.index.duplicated().any():
        joined = joined[~joined.index.duplicated(keep="first")]
    return (
        joined.groupby(["class", GRID_COL], observed=True)
        .size()
        .reset_index(name="n_obs")
    )


def _load_biel_grid_counts() -> pd.DataFrame:
    """Biel Mallorca flora → PAC/PAU cell counts (cuadrículas = PAC)."""
    if not BIEL_QUAD.exists() or not BIEL_FLORA.exists():
        print("Biel Excel missing — skip Biel layer")
        return pd.DataFrame(columns=["class", GRID_COL, "n_obs"])

    quad = pd.read_excel(BIEL_QUAD)
    if GRID_COL not in quad.columns:
        quad[GRID_COL] = quad["UTM"].map(utm_quad_to_code)
    pac_codes = set(
        quad.loc[quad["Illa"].eq("Mallorca") & quad[GRID_COL].notna(), GRID_COL]
        .astype(str)
        .unique()
    )

    flora = pd.read_excel(BIEL_FLORA)
    flora = flora.loc[flora["Illa"].eq("Mallorca")].copy()
    flora[GRID_COL] = flora["UTM 1x1"].map(parse_flora_utm_code)
    flora = flora.loc[flora[GRID_COL].notna()].copy()
    flora["class"] = np.where(flora[GRID_COL].isin(pac_codes), "PAC", "PAU")
    print(
        f"Biel Mallorca: PAC={int((flora['class']=='PAC').sum()):,} "
        f"PAU={int((flora['class']=='PAU').sum()):,} "
        f"pac_cells={len(pac_codes)}"
    )
    return (
        flora.groupby(["class", GRID_COL], observed=True)
        .size()
        .reset_index(name="n_obs")
    )


def _merge_counts(base: pd.DataFrame, biel: pd.DataFrame) -> pd.DataFrame:
    if biel.empty:
        return base
    both = pd.concat([base, biel], ignore_index=True)
    return (
        both.groupby(["class", GRID_COL], observed=True)["n_obs"]
        .sum()
        .reset_index()
    )


def main():
    apply_plot_style()
    net_all = gpd.read_file(FLORA_NET_PATH)[[GRID_COL, "geometry"]]
    net = _clip_mallorca(net_all)
    mallorca_codes = set(net[GRID_COL].astype(str))

    print(f"Loading {FULL_DATA_GPKG} …")
    counts = _merge_counts(_load_gpkg_grid_counts(net), _load_biel_grid_counts())
    counts = counts[counts[GRID_COL].astype(str).isin(mallorca_codes)].copy()
    for src in SOURCES:
        sub = counts[counts["class"] == src]
        print(f"  {src}: {int(sub['n_obs'].sum()):,} obs in {len(sub):,} cells")

    vmax = float(counts["n_obs"].quantile(0.99))
    print(f"Shared color scale vmax (p99): {vmax:.0f}")

    xlim = _mallorca_limits(net, counts[GRID_COL].unique())
    norm = Normalize(vmin=0, vmax=vmax)
    sm = ScalarMappable(cmap="YlOrRd", norm=norm)
    sm.set_array([])

    fig, axes = plt.subplots(2, 2, figsize=(11, 10), facecolor="white")
    for ax, src, letter in zip(axes.flat, SOURCES, PANEL_LETTERS):
        ax.set_facecolor("white")
        sub = counts[counts["class"] == src]
        m = net.merge(sub, on=GRID_COL, how="left")
        m["n_obs"] = m["n_obs"].fillna(0)

        m.plot(
            column="n_obs",
            ax=ax,
            cmap="YlOrRd",
            vmin=0,
            vmax=vmax,
            linewidth=0.05,
            edgecolor="#aaaaaa",
            legend=False,
        )
        ax.set_xlim(xlim[0], xlim[2])
        ax.set_ylim(xlim[1], xlim[3])
        ax.set_aspect("equal", adjustable="box")
        n_cells = len(sub)
        n_obs = int(sub["n_obs"].sum())
        # Panel letter only — narrative in manuscript caption
        ax.text(
            0.02,
            0.98,
            letter,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=14,
            fontweight="bold",
        )
        ax.text(
            0.98,
            0.02,
            f"{src}\nn={n_obs:,} · {n_cells:,} cells",
            transform=ax.transAxes,
            va="bottom",
            ha="right",
            fontsize=8,
        )
        ax.set_axis_off()

    fig.tight_layout(rect=[0, 0, 0.92, 1.0])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.75, pad=0.02)
    cbar.set_label("Observations per cell")

    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUTPUT_PDF, bbox_inches="tight", facecolor="white")
    print(f"Saved {OUTPUT_PNG}")
    print(f"Saved {OUTPUT_PDF}")
    plt.close(fig)


if __name__ == "__main__":
    main()
