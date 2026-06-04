"""
Train on mixes of PAU / POV / POU (0%, 50%, 100% of each source) and plot mean PAC AUC.
"""
from __future__ import annotations

import os
from itertools import product
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, TensorDataset, random_split

from models import MLP, BalancedBCELoss, DeepMaxEntLoss
from plot_style import apply as apply_plot_style
from run_model import _load_pac_eval, _load_source_subset

TRAIN_SOURCES = ["PAU", "POV", "POU"]
FRACTIONS = (0.0, 0.5, 1.0)
OUTPUT_CSV = Path("output/optimal_model_auc.csv")
OUTPUT_PNG = Path("output/optimal_model_auc.png")
SEED = 42
NUM_EPOCHS = 5


def _subsample(x: pd.DataFrame, y: pd.DataFrame, frac: float, seed: int):
    if frac <= 0 or len(x) == 0:
        return x.iloc[:0], y.iloc[:0]
    if frac >= 1.0:
        return x, y
    n = max(1, int(round(len(x) * frac)))
    idx = x.sample(n=min(n, len(x)), random_state=seed).index
    return x.loc[idx], y.loc[idx]


def _mix_label(pcts: dict[str, float]) -> str:
    return "_".join(f"{s}{int(pcts[s] * 100)}" for s in TRAIN_SOURCES)


def _build_training_mix(
    pcts: dict[str, float],
    species_cols: pd.Index,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame] | tuple[None, None]:
    xs, ys = [], []
    for i, src in enumerate(TRAIN_SOURCES):
        frac = pcts[src]
        if frac <= 0:
            continue
        x, y, _ = _load_source_subset(src, n_grids=None, seed=seed, grouped=False)
        x, y = _subsample(x, y, frac, seed=seed + i)
        if len(x) == 0:
            continue
        y = y.reindex(columns=species_cols, fill_value=0)
        xs.append(x)
        ys.append(y)

    if not xs:
        return None, None
    x_all = pd.concat(xs, axis=0)
    y_all = pd.concat(ys, axis=0)
    return x_all, y_all


def _mean_auc_on_pac(
    x_train: pd.DataFrame,
    y_train: pd.DataFrame,
    x_pac: pd.DataFrame,
    y_pac: pd.DataFrame,
    use_bce: bool,
) -> float:
    X_tr = x_train.values.astype(np.float32)
    Y_tr = y_train.values.astype(np.float32)
    X_pac = x_pac.values.astype(np.float32)

    scaler = StandardScaler().fit(X_tr)
    X_tr = scaler.transform(X_tr).astype(np.float32)
    X_pac = scaler.transform(X_pac).astype(np.float32)

    dataset = TensorDataset(torch.tensor(X_tr), torch.tensor(Y_tr))
    n = len(dataset)
    if n < 2:
        train_ds = val_ds = dataset
    else:
        n_train = max(1, int(0.8 * n))
        n_val = n - n_train
        if n_val == 0:
            n_train, n_val = n - 1, 1
        train_ds, val_ds = random_split(
            dataset,
            [n_train, n_val],
            generator=torch.Generator().manual_seed(SEED),
        )

    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=64)

    model = MLP(
        input_size=X_tr.shape[1],
        hidden_size=200,
        output_size=Y_tr.shape[1],
        hidden_layers=2,
    )
    criterion = BalancedBCELoss() if use_bce else DeepMaxEntLoss()
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=5, factor=0.5
    )

    best_val = float("inf")
    best_state = None

    def epoch(loader, train=True):
        model.train(train)
        total, count = 0.0, 0
        with torch.set_grad_enabled(train):
            for xb, yb in loader:
                loss = criterion(model(xb), yb)
                if train:
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                total += loss.item() * len(xb)
                count += len(xb)
        return total / max(count, 1)

    for _ in range(NUM_EPOCHS):
        epoch(train_loader, True)
        val_loss = epoch(val_loader, False)
        scheduler.step(val_loss)
        if val_loss < best_val:
            best_val = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

    if best_state is None:
        return float("nan")

    model.load_state_dict(best_state)
    model.eval()
    species = list(y_train.columns)
    with torch.no_grad():
        scores = model(torch.tensor(X_pac)).sigmoid().numpy()

    aucs = []
    for j, sp in enumerate(species):
        if sp not in y_pac.columns:
            continue
        try:
            aucs.append(roc_auc_score(y_pac[sp].values, scores[:, j]))
        except ValueError:
            pass
    return float(np.mean(aucs)) if aucs else float("nan")


def plot_optimal(results: pd.DataFrame, output_png: Path = OUTPUT_PNG) -> None:
    apply_plot_style()
    if results.empty:
        print("No results to plot.")
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    plot_df = results.sort_values("mean_auc_pac", ascending=True)
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(plot_df)))
    ax.barh(plot_df["label"], plot_df["mean_auc_pac"], color=colors)
    ax.set_xlabel("Mean per-species AUC on PAC")
    ax.set_title("Mixed training sources (fraction of rows per source)")
    ax.set_xlim(0, 1)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=150, bbox_inches="tight", facecolor="white")
    plt.show()
    print(f"Saved {output_png}")


def main(plot_only: bool = False):
    os.makedirs("output", exist_ok=True)

    if plot_only:
        if not OUTPUT_CSV.is_file():
            raise FileNotFoundError(
                f"No results at {OUTPUT_CSV}. Run without --plot-only first."
            )
        results = pd.read_csv(OUTPUT_CSV)
        print(f"Loaded {OUTPUT_CSV} ({len(results)} rows)")
        plot_optimal(results)
        return

    apply_plot_style()

    # Species shared across sources (intersection after full load)
    cols = None
    for src in TRAIN_SOURCES:
        _, y, _ = _load_source_subset(src, None, seed=SEED, grouped=False)
        cols = y.columns if cols is None else cols.intersection(y.columns)
    species_cols = pd.Index(cols)
    print(f"Training species (intersection): {len(species_cols)}")

    x_pac, y_pac = _load_pac_eval()
    rows = []

    combos = list(product(FRACTIONS, repeat=len(TRAIN_SOURCES)))
    for combo in combos:
        pcts = dict(zip(TRAIN_SOURCES, combo))
        if all(v <= 0 for v in pcts.values()):
            continue

        label = _mix_label(pcts)
        print(f"\n=== {label} ===")

        mixed = _build_training_mix(pcts, species_cols, seed=SEED)
        if mixed[0] is None:
            rows.append({**pcts, "label": label, "n_rows": 0, "mean_auc_pac": np.nan})
            continue

        x_train, y_train = mixed
        use_bce = pcts["PAU"] > 0
        mean_auc = _mean_auc_on_pac(x_train, y_train, x_pac, y_pac, use_bce=use_bce)
        print(f"  rows={len(x_train):,}  mean_auc={mean_auc:.4f}  loss={'BCE' if use_bce else 'MaxEnt'}")

        rows.append(
            {
                "label": label,
                "pct_PAU": pcts["PAU"],
                "pct_POV": pcts["POV"],
                "pct_POU": pcts["POU"],
                "n_rows": len(x_train),
                "mean_auc_pac": mean_auc,
            }
        )

    results = pd.DataFrame(rows).sort_values("mean_auc_pac", ascending=False)
    results.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved {OUTPUT_CSV}")
    plot_optimal(results)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Optimal mix of PAU/POV/POU sources.")
    p.add_argument(
        "--plot-only",
        action="store_true",
        help="Regenerate bar chart from output/optimal_model_auc.csv (skip training).",
    )
    main(plot_only=p.parse_args().plot_only)
