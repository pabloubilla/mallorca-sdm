#!/usr/bin/env python3
"""Q3 inventory: species accumulation under factorial mixes of PAU/POV/POU.

Subsample 0/50/100% of grid cells (PAU) or sites/cells (POV/POU) and count
unique species observed, plus how many PAC-rare species are recovered.

Uses MALLORCA_PROCESSED_DIR (default data/processed/curated).
PO compared in grouped (cell) mode for fair mix with PAU.
"""

from __future__ import annotations

import os
import sys
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plot_style import apply as apply_plot_style  # noqa: E402
from preprocessing import (  # noqa: E402
    PROCESSED_DIR,
    covariates_path,
    set_processed_dir,
    species_matrix_path,
)

OUT_DIR = ROOT / "output" / "curated"
TRAIN_SOURCES = ["PAU", "POV", "POU"]
FRACTIONS = (0.0, 0.5, 1.0)
SEED = 42


def _cell_matrix(src: str) -> pd.DataFrame:
    """One row per 1×1 km cell, binary species."""
    if src == "PAU":
        y = pd.read_csv(species_matrix_path("PAU", grouped=False), index_col="site_id")
        cov = pd.read_csv(covariates_path("PAU", grouped=False), index_col="site_id")
        if "UTMCODE1X1" in cov.columns:
            return y.groupby(cov["UTMCODE1X1"].astype(str)).max()
        return y
    return pd.read_csv(species_matrix_path(src, grouped=True), index_col="site_id")


def _pac_rare_species(processed_dir: str) -> set[str]:
    y = pd.read_csv(
        os.path.join(processed_dir, "species_matrix_PAC.csv"), index_col="site_id"
    )
    cov = pd.read_csv(
        os.path.join(processed_dir, "covariates_PAC.csv"), index_col="site_id"
    )
    if "UTMCODE1X1" in cov.columns:
        cell = y.groupby(cov["UTMCODE1X1"].astype(str)).max()
    else:
        cell = y
    pct = cell.sum(axis=0) / len(cell) * 100.0
    return set(pct[pct < 5.0].index)


def _subsample_rows(y: pd.DataFrame, frac: float, seed: int) -> pd.DataFrame:
    if frac <= 0 or len(y) == 0:
        return y.iloc[:0]
    if frac >= 1.0:
        return y
    n = max(1, int(round(len(y) * frac)))
    return y.sample(n=min(n, len(y)), random_state=seed)


def _species_in(y: pd.DataFrame) -> set[str]:
    if y.empty:
        return set()
    present = y.columns[y.sum(axis=0) > 0]
    return set(present)


def main() -> int:
    processed = os.environ.get("MALLORCA_PROCESSED_DIR", "data/processed/curated")
    set_processed_dir(processed)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_plot_style()

    mats = {src: _cell_matrix(src) for src in TRAIN_SOURCES}
    for src, y in mats.items():
        print(f"{src}: {len(y)} cells × {y.shape[1]} spp cols")
    rare = _pac_rare_species(PROCESSED_DIR)
    print(f"PAC-rare species (<5% cells): {len(rare)}")

    # Full single-source baselines + factorial mixes
    rows = []
    combos = list(product(FRACTIONS, repeat=len(TRAIN_SOURCES)))
    for i, combo in enumerate(combos):
        pcts = dict(zip(TRAIN_SOURCES, combo))
        if all(v <= 0 for v in pcts.values()):
            continue
        observed: set[str] = set()
        n_cells = 0
        for j, src in enumerate(TRAIN_SOURCES):
            frac = pcts[src]
            if frac <= 0:
                continue
            sub = _subsample_rows(mats[src], frac, seed=SEED + i * 10 + j)
            observed |= _species_in(sub)
            n_cells += len(sub)
        label = "_".join(f"{s}{int(pcts[s] * 100)}" for s in TRAIN_SOURCES)
        rows.append(
            {
                "label": label,
                "pct_PAU": pcts["PAU"],
                "pct_POV": pcts["POV"],
                "pct_POU": pcts["POU"],
                "n_cells": n_cells,
                "n_species": len(observed),
                "n_pac_rare_recovered": len(observed & rare),
                "pct_pac_rare_recovered": 100.0 * len(observed & rare) / len(rare)
                if rare
                else np.nan,
            }
        )

    # Also: exclusive contribution — species only in PAU vs only in PO
    pau_sp = _species_in(mats["PAU"])
    pov_sp = _species_in(mats["POV"])
    pou_sp = _species_in(mats["POU"])
    po_sp = pov_sp | pou_sp
    exclusives = {
        "only_PAU": len(pau_sp - po_sp),
        "only_PO": len(po_sp - pau_sp),
        "PAU_and_PO": len(pau_sp & po_sp),
        "only_POV_not_POU": len(pov_sp - pou_sp - pau_sp),
        "only_POU_not_POV": len(pou_sp - pov_sp - pau_sp),
    }
    pd.Series(exclusives).to_csv(OUT_DIR / "species_source_exclusives.csv", header=["n"])
    print("Exclusives:", exclusives)

    results = pd.DataFrame(rows).sort_values("n_species", ascending=False)
    results.to_csv(OUT_DIR / "species_accumulation_mix.csv", index=False)
    print(f"Wrote {OUT_DIR / 'species_accumulation_mix.csv'}")
    print(results.head(8).to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    plot_df = results.sort_values("n_species", ascending=True)
    axes[0].barh(plot_df["label"], plot_df["n_species"], color="#2c7fb8")
    axes[0].set_xlabel("Unique species observed")
    axes[0].set_title("Inventory mix (grouped cells)")
    axes[1].barh(plot_df["label"], plot_df["n_pac_rare_recovered"], color="#41ab5d")
    axes[1].set_xlabel("PAC-rare species recovered")
    axes[1].set_title("Rare recovery by mix")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "species_accumulation_mix.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
