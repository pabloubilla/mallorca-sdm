#!/usr/bin/env python3
"""A/B: PAU grouped MLP ± OpenTree Laplacian reg on species head → PAC AUC.

Reads laplacian from data/processed/curated/opentree/laplacian.npz
Writes output/curated/phylo_laplace_ab.csv
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, TensorDataset, random_split

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from models import MLP, BalancedBCELoss, phylo_laplacian_reg  # noqa: E402
from preprocessing import set_processed_dir  # noqa: E402
from run_model import (  # noqa: E402
    _aggregate_pixel_scores_to_cells,
    _load_pac_eval,
    _load_source_subset,
    _per_species_auc,
    _run_training_epoch,
)


def _align_L(L: np.ndarray, L_species: list[str], train_species: list[str]) -> np.ndarray:
    idx = {s: i for i, s in enumerate(L_species)}
    k = len(train_species)
    out = np.zeros((k, k), dtype=np.float32)
    mapped = []
    missing = 0
    for s in train_species:
        if s in idx:
            mapped.append(idx[s])
        else:
            mapped.append(None)
            missing += 1
    for i, mi in enumerate(mapped):
        if mi is None:
            continue
        for j, mj in enumerate(mapped):
            if mj is None:
                continue
            out[i, j] = L[mi, mj]
    if missing:
        print(f"  L align: {missing}/{k} train spp absent from Laplacian cache")
    return out


def train_eval_lambda(
    *,
    x_train: pd.DataFrame,
    y_train: pd.DataFrame,
    x_pac: pd.DataFrame,
    y_pac: pd.DataFrame,
    pac_grids: pd.Series | None,
    L_np: np.ndarray,
    lam: float,
    num_epochs: int,
    seed: int,
    out_dir: Path,
) -> dict:
    train_species = list(y_train.columns)
    X_train_bio = x_train.values.astype(np.float32)
    Y_train = y_train.values.astype(np.float32)
    X_pac_bio = x_pac.values.astype(np.float32)

    scaler = StandardScaler().fit(X_train_bio)
    X_train = scaler.transform(X_train_bio).astype(np.float32)
    X_pac = scaler.transform(X_pac_bio).astype(np.float32)

    X_train_t = torch.tensor(X_train)
    Y_train_t = torch.tensor(Y_train)
    L_t = torch.tensor(L_np, dtype=torch.float32)

    train_dataset = TensorDataset(X_train_t, Y_train_t)
    n = len(train_dataset)
    g = torch.Generator().manual_seed(seed)
    if n < 2:
        train_ds = val_ds = train_dataset
    else:
        n_train = max(1, int(0.8 * n))
        n_val = n - n_train
        if n_val == 0:
            n_train, n_val = n - 1, 1
        train_ds, val_ds = random_split(train_dataset, [n_train, n_val], generator=g)

    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=64)

    torch.manual_seed(seed)
    model = MLP(
        input_size=X_train.shape[1],
        hidden_size=200,
        output_size=Y_train.shape[1],
        hidden_layers=2,
    )
    criterion = BalancedBCELoss()
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=5, factor=0.5
    )

    ckpt = out_dir / f"best_model_PAU_grouped_phylo_lam{lam:g}.pt"
    best_val = float("inf")

    for epoch in range(1, num_epochs + 1):
        model.train()
        total_loss, total = 0.0, 0
        for xb, yb in train_loader:
            logits = model(xb)
            loss = criterion(logits, yb)
            if lam > 0:
                loss = loss + lam * phylo_laplacian_reg(model.output_layer.weight, L_t)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * len(xb)
            total += len(xb)
        train_loss = total_loss / max(total, 1)

        # val: BCE only (fair early-stop; phylo is train regularizer)
        val_loss = _run_training_epoch(
            model, criterion, optimizer, val_loader, train=False
        )
        scheduler.step(val_loss)
        if val_loss < best_val:
            best_val = val_loss
            torch.save(model.state_dict(), ckpt)
        print(
            f"  λ={lam:g} epoch {epoch}/{num_epochs}  "
            f"train={train_loss:.4f} val={val_loss:.4f}"
        )

    model.load_state_dict(torch.load(ckpt, weights_only=True))
    model.eval()
    with torch.no_grad():
        y_score = model(torch.tensor(X_pac)).sigmoid().numpy()

    y_score_df = pd.DataFrame(y_score, columns=train_species, index=y_pac.index)
    if pac_grids is not None:
        y_pac_eval, y_score_df = _aggregate_pixel_scores_to_cells(
            y_pac, y_score_df, pac_grids
        )
    else:
        y_pac_eval = y_pac

    # AUC only on species present in both train head and PAC labels
    eval_spp = [s for s in train_species if s in y_pac_eval.columns]
    auc_scores = _per_species_auc(y_pac_eval, y_score_df, eval_spp)
    auc = pd.Series(auc_scores)
    valid = auc.dropna()
    return {
        "lambda": lam,
        "mean_auc": float(valid.mean()) if len(valid) else float("nan"),
        "median_auc": float(valid.median()) if len(valid) else float("nan"),
        "n_valid_auc": int(valid.notna().sum()),
        "n_species_train": len(train_species),
        "best_val_loss": float(best_val),
        "n_train_rows": len(y_train),
        "auc_series": auc,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-dir", default="data/processed/curated")
    parser.add_argument(
        "--laplace",
        default=None,
        help="Path to laplacian.npz (default <processed-dir>/opentree/laplacian.npz)",
    )
    parser.add_argument(
        "--lambdas",
        default="0,1e-3,1e-2,1e-1",
        help="Comma-separated λ values",
    )
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-dir", default="output/curated")
    args = parser.parse_args()

    set_processed_dir(args.processed_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lap_path = Path(
        args.laplace or os.path.join(args.processed_dir, "opentree", "laplacian.npz")
    )
    if not lap_path.is_file():
        raise SystemExit(
            f"Missing {lap_path}. Run scripts/build_opentree_laplacian.py first."
        )

    z = np.load(lap_path, allow_pickle=True)
    L_full = z["L"]
    L_species = [str(s) for s in z["species"].tolist()]
    sigma = float(z["sigma"])
    n_matched = int(z["n_matched"])
    match_rate = float(z["match_rate"]) if "match_rate" in z.files else float("nan")
    print(
        f"Laplacian {lap_path}: K={L_full.shape[0]} sigma={sigma:.4g} "
        f"n_matched={n_matched} match_rate={match_rate:.1%}"
    )

    print("Loading PAU grouped…")
    x_train, y_train, meta = _load_source_subset(
        "PAU", n_grids=None, seed=args.seed, grouped=True
    )
    x_pac, y_pac, pac_grids = _load_pac_eval()
    # Align PAC columns to train species (reindex fill 0 for missing)
    y_pac = y_pac.reindex(columns=list(y_train.columns), fill_value=0)

    L_np = _align_L(L_full, L_species, list(y_train.columns))
    print(
        f"Train cells={meta['n_grids']} rows={len(y_train)} spp={y_train.shape[1]} "
        f"| PAC pixels={len(y_pac)}"
    )

    lambdas = [float(x) for x in args.lambdas.split(",") if x.strip()]
    rows = []
    auc_by_lam: dict[float, pd.Series] = {}
    for lam in lambdas:
        print(f"\n=== λ = {lam:g} ===")
        res = train_eval_lambda(
            x_train=x_train,
            y_train=y_train,
            x_pac=x_pac,
            y_pac=y_pac,
            pac_grids=pac_grids,
            L_np=L_np,
            lam=lam,
            num_epochs=args.epochs,
            seed=args.seed,
            out_dir=out_dir,
        )
        auc_by_lam[lam] = res.pop("auc_series")
        rows.append(
            {
                **{k: v for k, v in res.items()},
                "sigma": sigma,
                "n_matched_tips": n_matched,
                "match_rate": match_rate,
            }
        )
        print(
            f"  → mean AUC {res['mean_auc']:.4f}  median {res['median_auc']:.4f}  "
            f"n_valid={res['n_valid_auc']}"
        )

    summary = pd.DataFrame(rows)
    # Same species validity: intersect non-NaN across all λ for delta fairness
    common_valid = None
    for s in auc_by_lam.values():
        mask = s.notna()
        common_valid = mask if common_valid is None else (common_valid & mask)
    if common_valid is not None and common_valid.any():
        control = auc_by_lam[lambdas[0]]
        for i, lam in enumerate(lambdas):
            a = auc_by_lam[lam][common_valid]
            summary.loc[i, "mean_auc_common"] = float(a.mean())
            summary.loc[i, "median_auc_common"] = float(a.median())
            summary.loc[i, "delta_mean_vs_control"] = float(
                a.mean() - control[common_valid].mean()
            )
        summary["n_common_valid"] = int(common_valid.sum())

    out_csv = out_dir / "phylo_laplace_ab.csv"
    summary.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")
    print(summary.to_string(index=False))

    control_mean = summary.loc[summary["lambda"] == 0.0, "mean_auc_common"]
    if len(control_mean):
        c = float(control_mean.iloc[0])
        best = summary.loc[summary["lambda"] > 0].sort_values(
            "mean_auc_common", ascending=False
        )
        if len(best):
            b = best.iloc[0]
            ok = b["mean_auc_common"] >= c + 0.005
            print(
                f"\nVerdict: best λ={b['lambda']:g} mean_common={b['mean_auc_common']:.4f} "
                f"vs control={c:.4f} (Δ={b['delta_mean_vs_control']:+.4f}) "
                f"{'PASS (+≥0.005)' if ok else 'NO GAIN'}"
            )


if __name__ == "__main__":
    main()
