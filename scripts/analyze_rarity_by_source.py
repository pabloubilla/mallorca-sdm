#!/usr/bin/env python3
"""Q1 descriptive: rarity by source; does PO under-sample PAC-common species?

Uses PAC prevalence bins (<5%, 5–30%, >30%) matching sampling_effort.py.
Reads matrices under MALLORCA_PROCESSED_DIR (default data/processed/curated).

Outputs under output/curated/:
  rarity_by_source.csv
  rarity_by_source.png
  common_pac_in_po_rarification.csv
  common_pac_in_po_rarification.png
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plot_style import apply as apply_plot_style  # noqa: E402
from preprocessing import (  # noqa: E402
    PROCESSED_DIR,
    species_matrix_path,
    set_processed_dir,
)

OUT_DIR = ROOT / "output" / "curated"
TRAIN_SOURCES = ("PAU", "POV", "POU")
BINS = (
    ("rare", 0.0, 5.0),
    ("medium", 5.0, 30.0),
    ("common", 30.0, 100.1),
)


def _pac_prevalence(processed_dir: str) -> pd.Series:
    y = pd.read_csv(
        os.path.join(processed_dir, "species_matrix_PAC.csv"), index_col="site_id"
    )
    # Aggregate pixel rows to cells if covariates carry UTMCODE1X1
    cov = pd.read_csv(
        os.path.join(processed_dir, "covariates_PAC.csv"), index_col="site_id"
    )
    if "UTMCODE1X1" in cov.columns:
        grids = cov["UTMCODE1X1"].astype(str)
        cell = y.groupby(grids).max()
    else:
        cell = y
    n = len(cell)
    return (cell.sum(axis=0) / n * 100.0).rename("pac_prevalence_pct")


def _bin_label(pct: float) -> str:
    for name, lo, hi in BINS:
        if lo <= pct < hi:
            return name
    return "common"


def _source_matrix(src: str, grouped: bool) -> pd.DataFrame:
    path = species_matrix_path(src, grouped=grouped)
    return pd.read_csv(path, index_col="site_id")


def main() -> int:
    processed = os.environ.get("MALLORCA_PROCESSED_DIR", "data/processed/curated")
    set_processed_dir(processed)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_plot_style()

    prev = _pac_prevalence(PROCESSED_DIR)
    prev_df = prev.reset_index()
    prev_df.columns = ["scientificName", "pac_prevalence_pct"]
    prev_df["pac_bin"] = prev_df["pac_prevalence_pct"].map(_bin_label)
    print(prev_df["pac_bin"].value_counts().to_string())

    rows = []
    for src in TRAIN_SOURCES:
        for grouped, mode in ((False, "point"), (True, "grouped")):
            if src == "PAU" and grouped:
                continue  # PA already grid
            y = _source_matrix(src, grouped=(src != "PAU" and grouped))
            if src == "PAU":
                cov = pd.read_csv(
                    os.path.join(PROCESSED_DIR, "covariates_PAU.csv"),
                    index_col="site_id",
                )
                if "UTMCODE1X1" in cov.columns:
                    y = y.groupby(cov["UTMCODE1X1"].astype(str)).max()
                mode = "grid"

            present = y.columns[y.sum(axis=0) > 0]
            occ = y.sum(axis=0)
            for sp in present:
                if sp not in prev.index:
                    continue
                rows.append(
                    {
                        "source": src,
                        "mode": mode,
                        "scientificName": sp,
                        "n_sites_present": int(occ[sp]),
                        "pac_prevalence_pct": float(prev[sp]),
                        "pac_bin": _bin_label(float(prev[sp])),
                    }
                )

    sp_src = pd.DataFrame(rows)
    sp_src.to_csv(OUT_DIR / "rarity_species_by_source.csv", index=False)

    # Share of species and of occurrence-sites by PAC bin
    summary = []
    for (src, mode), sub in sp_src.groupby(["source", "mode"]):
        n_sp = len(sub)
        w = sub["n_sites_present"]
        for bin_name, _, _ in BINS:
            b = sub[sub["pac_bin"] == bin_name]
            summary.append(
                {
                    "source": src,
                    "mode": mode,
                    "pac_bin": bin_name,
                    "n_species": len(b),
                    "pct_species": 100.0 * len(b) / n_sp if n_sp else np.nan,
                    "n_site_presences": int(b["n_sites_present"].sum()),
                    "pct_site_presences": 100.0 * b["n_sites_present"].sum() / w.sum()
                    if w.sum()
                    else np.nan,
                }
            )
    # PAC reference
    for bin_name, _, _ in BINS:
        b = prev_df[prev_df["pac_bin"] == bin_name]
        summary.append(
            {
                "source": "PAC",
                "mode": "grid",
                "pac_bin": bin_name,
                "n_species": len(b),
                "pct_species": 100.0 * len(b) / len(prev_df),
                "n_site_presences": np.nan,
                "pct_site_presences": np.nan,
            }
        )
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(OUT_DIR / "rarity_by_source.csv", index=False)
    print(f"Wrote {OUT_DIR / 'rarity_by_source.csv'}")

    # Plot species % by bin
    plot_df = summary_df[summary_df["mode"].isin(["point", "grid"])].copy()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    sources = ["PAC", "PAU", "POV", "POU"]
    x = np.arange(len(sources))
    width = 0.25
    for i, bin_name in enumerate(["rare", "medium", "common"]):
        vals = []
        for src in sources:
            row = plot_df[
                (plot_df["source"] == src)
                & (plot_df["pac_bin"] == bin_name)
                & (
                    (plot_df["mode"] == "grid")
                    if src in ("PAC", "PAU")
                    else (plot_df["mode"] == "point")
                )
            ]
            vals.append(float(row["pct_species"].iloc[0]) if len(row) else 0.0)
        ax.bar(x + (i - 1) * width, vals, width, label=bin_name)
    ax.set_xticks(x)
    ax.set_xticklabels(sources)
    ax.set_ylabel("% of species (by PAC prevalence bin)")
    ax.set_title("Species composition by PAC rarity bin and source (curated)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT_DIR / "rarity_by_source.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Q1b: PAC-common species — rarification in PO vs PAU at cell grain.
    # PO must be grouped: point-site % is not comparable to PAC cell %.
    src_mats: dict[str, pd.DataFrame] = {}
    for src in TRAIN_SOURCES:
        if src == "PAU":
            y = _source_matrix(src, grouped=False)
            cov = pd.read_csv(
                os.path.join(PROCESSED_DIR, "covariates_PAU.csv"),
                index_col="site_id",
            )
            if "UTMCODE1X1" in cov.columns:
                y = y.groupby(cov["UTMCODE1X1"].astype(str)).max()
        else:
            y = _source_matrix(src, grouped=True)
        src_mats[src] = y
        print(f"loaded {src} (cells): {y.shape}")

    common = prev_df[prev_df["pac_bin"] == "common"]["scientificName"]
    rar_rows = []
    for sp in common:
        pac_pct = float(prev[sp])
        row = {"scientificName": sp, "pac_prevalence_pct": pac_pct}
        for src in TRAIN_SOURCES:
            y = src_mats[src]
            if sp not in y.columns:
                row[f"{src}_present"] = False
                row[f"{src}_pct_sites"] = 0.0
                row[f"{src}_rarify_ratio"] = np.inf
                continue
            n_sites = len(y)
            pct = 100.0 * (y[sp] > 0).sum() / n_sites if n_sites else 0.0
            row[f"{src}_present"] = bool((y[sp] > 0).any())
            row[f"{src}_pct_sites"] = pct
            # rarification ratio: PAC prevalence / source prevalence ( >1 = looks rarer in source)
            row[f"{src}_rarify_ratio"] = pac_pct / pct if pct > 0 else np.inf
        rar_rows.append(row)
    rar_df = pd.DataFrame(rar_rows)
    rar_df.to_csv(OUT_DIR / "common_pac_in_po_rarification.csv", index=False)

    # Summary: median rarify ratio among PAC-common spp that appear in source
    print("\nPAC-common species rarification (median PAC%/source% among present):")
    for src in TRAIN_SOURCES:
        col = f"{src}_rarify_ratio"
        finite = rar_df[np.isfinite(rar_df[col])][col]
        missing = (~rar_df[f"{src}_present"]).sum()
        print(
            f"  {src}: n_common={len(rar_df)} absent={missing} "
            f"median_ratio={finite.median():.2f} mean_ratio={finite.mean():.2f}"
        )

    fig, ax = plt.subplots(figsize=(7, 4.5))
    data, labels = [], []
    for src in TRAIN_SOURCES:
        col = f"{src}_rarify_ratio"
        vals = rar_df[np.isfinite(rar_df[col])][col].clip(upper=50)
        data.append(vals)
        labels.append(src)
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.set_ylabel("PAC% / source% (clipped ≤50)")
    ax.set_title("Do sources make PAC-common species look rarer? (ratio>1)")
    fig.tight_layout()
    fig.savefig(
        OUT_DIR / "common_pac_in_po_rarification.png", dpi=150, bbox_inches="tight"
    )
    plt.close(fig)
    print(f"Wrote {OUT_DIR / 'common_pac_in_po_rarification.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
