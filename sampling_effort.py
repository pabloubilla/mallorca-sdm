"""Learning curves by training size, stratified by species prevalence on PAC."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from preprocessing import grid_cells_path, is_po_source
from run_model import MIN_TRAIN_SITES, run_model

# Training units = COD10X10 grid cells (see preprocessing.py), not raw GPS points
SIZE_GRID = [1, 10,100,200,500,1000,1600]
TRAIN_SOURCES = ["PAU", "POU", "POV"]
EVAL_SOURCE = "PAC"

# "all" = todas las especies (comportamiento original)
# "intersection_pac_pov_pou" = especies con >=1 presencia en PAC, POV y POU
EVAL_SPECIES_MODE_DEFAULT = "all"
EVAL_SPECIES_MODES = ("all", "intersection_pac_pov_pou")

# Prevalence on PAC: % of COD10X10 cells with presence (not training n_grids).
PREVALENCE_BINS = [
    ("rare", "<5% PAC cells", lambda p: p < 5),
    ("medium", "5–30% PAC cells", lambda p: (p >= 5) & (p <= 30)),
    ("common", ">30% PAC cells", lambda p: p > 30),
]


def _output_paths(eval_species_mode: str) -> tuple[Path, Path]:
    suffix = "" if eval_species_mode == "all" else f"_{eval_species_mode}"
    return (
        Path(f"output/sampling_effort_auc{suffix}.csv"),
        Path(f"output/sampling_effort_auc_by_prevalence{suffix}.png"),
    )


def _species_with_presence(source: str) -> pd.Index:
    """Presence counts on the natural unit (points for PO, grids for PA)."""
    y = pd.read_csv(
        f"data/processed/species_matrix_{source}.csv", index_col="site_id"
    )
    counts = y.sum(axis=0)
    return counts.index[counts > 0]


def _n_grids_available(source: str) -> int:
    return len(pd.read_csv(grid_cells_path(source)))


def resolve_eval_species(eval_species_mode: str) -> pd.Index | None:
    """
    Especies sobre las que agregar AUC en evaluación.
    None = todas (sin filtro adicional).
    """
    if eval_species_mode == "all":
        return None
    if eval_species_mode == "intersection_pac_pov_pou":
        sets = [
            set(_species_with_presence(s)) for s in ("PAC", "POV", "POU")
        ]
        common = sets[0].intersection(*sets[1:])
        return pd.Index(sorted(common))
    raise ValueError(
        f"eval_species_mode must be one of {EVAL_SPECIES_MODES}, got {eval_species_mode!r}"
    )


def _pac_prevalence_pct() -> tuple[pd.Series, int]:
    """Per-species % of PAC grid cells with presence; returns (pct, n_grids)."""
    y = pd.read_csv(
        f"data/processed/species_matrix_{EVAL_SOURCE}.csv", index_col="site_id"
    )
    n_grids = len(y)
    if n_grids == 0:
        raise ValueError("PAC species matrix has no grid rows")
    counts = y.sum(axis=0)
    return counts / n_grids * 100.0, n_grids


def _species_in_bin(prevalence: pd.Series, mask_fn) -> pd.Index:
    return prevalence.index[mask_fn(prevalence.values)]


def _restrict_species(species: pd.Index, eval_species: pd.Index | None) -> pd.Index:
    if eval_species is None:
        return species
    return species.intersection(eval_species)


def _sizes_for_source(n_grids: int) -> list[int]:
    return [s for s in SIZE_GRID if MIN_TRAIN_SITES <= s <= n_grids]


def _mean_auc_in_bin(auc: pd.Series, species: pd.Index) -> float:
    sub = auc.reindex(species).dropna()
    return float(sub.mean()) if len(sub) else np.nan


def sampling_effort_per_source(
    eval_species_mode: str = EVAL_SPECIES_MODE_DEFAULT,
) -> pd.DataFrame:
    """
    Train at increasing sizes; plot mean PAC AUC vs. size in 3 panels by prevalence bin.

    eval_species_mode:
        - "all": evaluar en todas las especies (bins de prevalencia en PAC completos).
        - "intersection_pac_pov_pou": solo especies con presencia en PAC, POV y POU.
    """
    eval_species = resolve_eval_species(eval_species_mode)
    output_csv, output_png = _output_paths(eval_species_mode)

    prevalence_pct, n_pac_grids = _pac_prevalence_pct()
    bin_species = {
        bid: _restrict_species(_species_in_bin(prevalence_pct, fn), eval_species)
        for bid, _label, fn in PREVALENCE_BINS
    }

    print(f"Eval species mode: {eval_species_mode}")
    print(
        f"PAC prevalence: % of {n_pac_grids} COD10X10 cells "
        f"(max {prevalence_pct.max():.1f}%)"
    )
    if eval_species is not None:
        print(f"  species in evaluation set: {len(eval_species)}")
    for bid, label, _ in PREVALENCE_BINS:
        print(f"  Bin {label}: {len(bin_species[bid])} species")

    rows: list[dict] = []

    for source in TRAIN_SOURCES:
        n_max = _n_grids_available(source)
        sizes = _sizes_for_source(n_max)
        if not sizes:
            print(
                f"[{source}] skip: fewer than {MIN_TRAIN_SITES} grids ({n_max})"
            )
            continue

        unit_note = "grids → points inside" if is_po_source(source) else "grids"
        print(f"\n=== {source}: sizes {sizes} (max grids {n_max}, {unit_note}) ===")
        for size in sizes:
            print(f"--- {source}, n_grids={size} ---")
            auc_by_sp = run_model(source, size_train=size)
            if eval_species is not None:
                auc_by_sp = auc_by_sp.reindex(eval_species)

            for bin_id, bin_label, _ in PREVALENCE_BINS:
                sp = bin_species[bin_id]
                rows.append(
                    {
                        "eval_species_mode": eval_species_mode,
                        "train_source": source,
                        "n_grids": size,
                        "n_train": size,  # legacy column = grids (effort axis)
                        "prevalence_bin": bin_id,
                        "prevalence_label": bin_label,
                        "n_species_bin": len(sp),
                        "n_auc_computed": int(
                            auc_by_sp.reindex(sp).notna().sum()
                        ),
                        "mean_auc_pac": _mean_auc_in_bin(auc_by_sp, sp),
                    }
                )

    results = pd.DataFrame(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_csv, index=False)
    print(f"\nSaved {output_csv}")

    if results.empty:
        return results

    title_suffix = (
        ""
        if eval_species_mode == "all"
        else f" | eval: {eval_species_mode}"
    )
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    for ax, (bin_id, bin_label, _) in zip(axes, PREVALENCE_BINS):
        sub_bin = results[results["prevalence_bin"] == bin_id]
        for source in TRAIN_SOURCES:
            sub = sub_bin[sub_bin["train_source"] == source].sort_values("n_train")
            if sub.empty:
                continue
            ax.plot(
                sub["n_train"],
                sub["mean_auc_pac"],
                marker="o",
                label=source,
            )
        ax.set_xscale("log")
        ax.set_xlabel("Training grids COD10X10 (log scale)")
        ax.set_title(f"{EVAL_SOURCE} prevalence: {bin_label}")
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Mean per-species AUC")
    axes[-1].legend(loc="lower right")
    fig.suptitle(
        "Sampling effort vs. transfer performance by species prevalence"
        + title_suffix,
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved {output_png}")
    return results


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Sampling effort curves (PAC eval).")
    p.add_argument(
        "--eval-species",
        choices=EVAL_SPECIES_MODES,
        default=EVAL_SPECIES_MODE_DEFAULT,
        help=(
            "Species set for aggregating AUC: 'all' (default) or "
            "'intersection_pac_pov_pou' (>=1 presence in PAC, POV and POU)."
        ),
    )
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    sampling_effort_per_source(eval_species_mode=args.eval_species)
