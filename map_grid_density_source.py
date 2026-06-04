"""2×2 map: observation density per UTMCODE1X1 cell for PAC, PAU, POV, POU."""
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from plot_style import apply as apply_plot_style
from preprocessing import EXTRACTED_CSV, FLORA_NET_PATH, GRID_COL, assign_cod10_grid

SOURCES = ["PAC", "PAU", "POV", "POU"]
OUTPUT_PNG = Path("output/grid_density_by_source.png")

## increase the font size
plt.rcParams.update({'font.size': 14})


def _mallorca_limits(net: gpd.GeoDataFrame, active_codes: pd.Index) -> tuple[float, float, float, float]:
    """Extent of grid cells with data, padded (~Mallorca focus)."""
    sub = net[net[GRID_COL].isin(active_codes)]
    if sub.empty:
        return net.total_bounds
    xmin, ymin, xmax, ymax = sub.total_bounds
    pad_x = (xmax - xmin) * 0.05 or 5000
    pad_y = (ymax - ymin) * 0.05 or 5000
    return xmin - pad_x, ymin - pad_y, xmax + pad_x, ymax + pad_y


def main():
    apply_plot_style()
    net = gpd.read_file(FLORA_NET_PATH)[[GRID_COL, "geometry"]]

    df = pd.read_csv(EXTRACTED_CSV)
    df = df[df["slope"].notna()].copy()
    print(f"Records with slope: {len(df):,}")
    df = assign_cod10_grid(df)

    counts = df.groupby(["class", GRID_COL]).size().reset_index(name="n_obs")
    vmax = float(counts["n_obs"].quantile(0.99))
    print(f"Shared color scale vmax (p99): {vmax:.0f}")

    xlim = _mallorca_limits(net, counts[GRID_COL].unique())
    norm = Normalize(vmin=0, vmax=vmax)
    sm = ScalarMappable(cmap="YlOrRd", norm=norm)
    sm.set_array([])

    fig, axes = plt.subplots(2, 2, figsize=(11, 10), facecolor="white")
    for ax, src in zip(axes.flat, SOURCES):
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
        ax.set_title(f"{src} — {n_obs:,} obs in {n_cells:,} cells")
        ax.set_axis_off()

    fig.suptitle(f"", fontsize=14, y=0.98)
    fig.tight_layout(rect=[0, 0, 0.92, 0.96])
    cbar = fig.colorbar(sm, ax=axes.ravel().tolist(), shrink=0.75, pad=0.02)
    cbar.set_label("n° observations per cell (species in PA and observations in PO)")

    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=150, bbox_inches="tight", facecolor="white")
    print(f"Saved {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()
