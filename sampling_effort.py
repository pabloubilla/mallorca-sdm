"""Learning curves: point vs grouped PO; 3×3 plot (effort rows × prevalence cols)."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_style import apply as apply_plot_style
from preprocessing import grid_cells_path, species_matrix_path
from run_model import MIN_TRAIN_SITES, run_model

SIZE_GRID = [1, 10, 100, 200, 500, 1000, 1600]
TRAIN_SOURCES = ["PAU", "POU", "POV"]
EVAL_SOURCE = "PAC"
EVAL_SPECIES_MODE_DEFAULT = "all"
EVAL_SPECIES_MODES = ("all", "intersection_pac_pov_pou")

PREVALENCE_BINS = [
    ("rare", "<5% PAC cells", lambda p: p < 5),
    ("medium", "5–30% PAC cells", lambda p: (p >= 5) & (p <= 30)),
    ("common", ">30% PAC cells", lambda p: p > 30),
]

# (row, data_mode filter, x column, x label)
PLOT_ROWS = [
    ("point", "point", "n_grids", "1×1 km grids (point mode)"),
    ("point", "point", "n_train_rows", "Training rows (point mode)"),
    ("grouped", "grouped", "n_grids", "1×1 km grids (PO grouped)"),
]


def _output_paths(eval_species_mode: str) -> tuple[Path, Path]:
    suffix = "" if eval_species_mode == "all" else f"_{eval_species_mode}"
    return (
        Path(f"output/sampling_effort_auc{suffix}.csv"),
        Path(f"output/sampling_effort_auc_by_prevalence{suffix}.png"),
    )


def resolve_eval_species(eval_species_mode: str) -> pd.Index | None:
    if eval_species_mode == "all":
        return None
    if eval_species_mode == "intersection_pac_pov_pou":
        sets = []
        for s in ("PAC", "POV", "POU"):
            counts = pd.read_csv(
                species_matrix_path(s, grouped=False), index_col="site_id"
            ).sum(axis=0)
            sets.append(set(counts.index[counts > 0]))
        return pd.Index(sorted(sets[0].intersection(*sets[1:])))
    raise ValueError(f"Unknown eval_species_mode: {eval_species_mode!r}")


def plot_sampling_effort(
    results: pd.DataFrame,
    output_png: Path,
    eval_species_mode: str,
) -> None:
    apply_plot_style()
    if results.empty:
        print("No results to plot.")
        return

    title_suffix = (
        "" if eval_species_mode == "all" else f" | {eval_species_mode}"
    )
    fig, axes = plt.subplots(3, 3, figsize=(15, 11), sharey="row")
    for col, (bin_id, bin_label, _) in enumerate(PREVALENCE_BINS):
        sub_bin = results[results["prevalence_bin"] == bin_id]
        for row, (mode, _m, x_col, x_label) in enumerate(PLOT_ROWS):
            ax = axes[row, col]
            sub = sub_bin[sub_bin["data_mode"] == mode]
            for source in TRAIN_SOURCES:
                s = sub[sub["train_source"] == source].sort_values(x_col)
                if s.empty:
                    continue
                ax.plot(s[x_col], s["mean_auc_pac"], marker="o", label=source)
            ax.set_xscale("log")
            ax.set_xlabel(x_label, fontsize=8)
            if row == 0:
                ax.set_title(bin_label)
            ax.grid(True, alpha=0.3)
        axes[row, 0].set_ylabel("Mean AUC")

    axes[2, 2].legend(loc="lower right", fontsize=8)
    fig.suptitle("Sampling effort 3×3: grids / rows / PO grouped" + title_suffix)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved {output_png}")


def sampling_effort_per_source(
    eval_species_mode: str = EVAL_SPECIES_MODE_DEFAULT,
    plot_only: bool = False,
) -> pd.DataFrame:
    output_csv, output_png = _output_paths(eval_species_mode)

    if plot_only:
        if not output_csv.is_file():
            raise FileNotFoundError(
                f"No results at {output_csv}. Run without --plot-only first."
            )
        results = pd.read_csv(output_csv)
        print(f"Loaded {output_csv} ({len(results)} rows)")
        plot_sampling_effort(results, output_png, eval_species_mode)
        return results

    eval_species = resolve_eval_species(eval_species_mode)

    y_pac = pd.read_csv(
        f"data/processed/species_matrix_{EVAL_SOURCE}.csv", index_col="site_id"
    )
    n_pac_grids = len(y_pac)
    prevalence_pct = y_pac.sum(axis=0) / n_pac_grids * 100.0
    bin_species = {
        bid: (
            prevalence_pct.index[fn(prevalence_pct.values)]
            if eval_species is None
            else prevalence_pct.index[fn(prevalence_pct.values)].intersection(
                eval_species
            )
        )
        for bid, _label, fn in PREVALENCE_BINS
    }

    print(f"Eval species mode: {eval_species_mode}")
    print(f"PAC prevalence: % of {n_pac_grids} cells (max {prevalence_pct.max():.1f}%)")
    for bid, label, _ in PREVALENCE_BINS:
        print(f"  {label}: {len(bin_species[bid])} species")

    rows: list[dict] = []

    for grouped in (False, True):
        mode = "grouped" if grouped else "point"
        print(f"\n######## data_mode={mode} ########")

        for source in TRAIN_SOURCES:
            n_max = len(pd.read_csv(grid_cells_path(source, grouped=grouped)))
            sizes = [s for s in SIZE_GRID if MIN_TRAIN_SITES <= s <= n_max]
            if not sizes:
                print(f"[{source}] skip ({mode}): max grids {n_max}")
                continue

            print(f"=== {source} ({mode}) sizes {sizes} ===")
            for size in sizes:
                print(f"--- {source} {mode} n_grids={size} ---")
                auc_by_sp, n_train_rows = run_model(
                    source, size_train=size, grouped=grouped
                )
                if eval_species is not None:
                    auc_by_sp = auc_by_sp.reindex(eval_species)

                for bin_id, bin_label, _ in PREVALENCE_BINS:
                    sp = bin_species[bin_id]
                    rows.append(
                        {
                            "eval_species_mode": eval_species_mode,
                            "data_mode": mode,
                            "train_source": source,
                            "n_grids": size,
                            "n_train_rows": n_train_rows,
                            "prevalence_bin": bin_id,
                            "prevalence_label": bin_label,
                            "n_species_bin": len(sp),
                            "n_auc_computed": int(
                                auc_by_sp.reindex(sp).notna().sum()
                            ),
                            "mean_auc_pac": float(
                                auc_by_sp.reindex(sp).dropna().mean()
                            )
                            if auc_by_sp.reindex(sp).notna().any()
                            else np.nan,
                        }
                    )

    results = pd.DataFrame(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_csv, index=False)
    print(f"\nSaved {output_csv}")

    plot_sampling_effort(results, output_png, eval_species_mode)
    return results


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Sampling effort (3×3 plot).")
    p.add_argument(
        "--eval-species",
        choices=EVAL_SPECIES_MODES,
        default=EVAL_SPECIES_MODE_DEFAULT,
    )
    p.add_argument(
        "--plot-only",
        action="store_true",
        help="Regenerate 3×3 figure from existing output CSV (skip training).",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    sampling_effort_per_source(
        eval_species_mode=args.eval_species,
        plot_only=args.plot_only,
    )
