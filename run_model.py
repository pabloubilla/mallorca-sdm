from models import MLP, BalancedBCELoss, DeepMaxEntLoss
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import BCEWithLogitsLoss
from torch.utils.data import DataLoader, TensorDataset, random_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from typing import List, Dict
import os

from preprocessing import (
    GRID_COL,
    MIN_SPECIES_LIST_PAC_GRIDS,
    MIN_TRAIN_OCC_GRID,
    MIN_TRAIN_OCC_POINT,
    covariates_path,
    file_suffix,
    grid_cells_path,
    is_po_source,
    species_matrix_path,
)

MIN_TRAIN_SITES = 1  # min COD10X10 cells when subsampling training effort


def _filter_species_by_occurrence(
    y: pd.DataFrame, min_occ: int, unit_label: str
) -> pd.DataFrame:
    counts = y.sum(axis=0)
    keep = counts[counts >= min_occ].index
    dropped = y.shape[1] - len(keep)
    if dropped:
        print(
            f"  species filter ({unit_label}, >={min_occ}): "
            f"keep {len(keep)}, drop {dropped}"
        )
    return y.loc[:, keep]


def _sample_grids(rng_grids: pd.Index, n_grids: int, seed: int) -> pd.Index:
    if n_grids >= len(rng_grids):
        return rng_grids
    picked = pd.Series(rng_grids.astype(str).values).sample(
        n=n_grids, random_state=seed
    )
    return pd.Index(picked.values)


def _load_source_subset(
    source: str,
    n_grids: int | None,
    seed: int = 42,
    min_occ_grid: int = MIN_TRAIN_OCC_GRID,
    min_occ_point: int = MIN_TRAIN_OCC_POINT,
    grouped: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Load (X, Y) for training.

    grouped=False: PAU grids; POV/POU GPS points (subsample cells, keep points).
    grouped=True: all sources at grid cell (PO mean bioclim per cell).
    """
    meta: dict = {"n_grids": None, "n_points": None, "unit": "grid", "grouped": grouped}

    y = pd.read_csv(species_matrix_path(source, grouped), index_col="site_id")
    x = pd.read_csv(covariates_path(source, grouped), index_col="site_id")
    bioclim_cols = [c for c in x.columns if c != GRID_COL]
    common = y.index.intersection(x.index)
    y, x = y.loc[common], x.loc[common]

    use_points = is_po_source(source) and not grouped

    if use_points:
        meta["unit"] = "point"
        y = _filter_species_by_occurrence(y, min_occ_point, "points (full PO)")
        grids_available = pd.read_csv(grid_cells_path(source, grouped=False))[
            GRID_COL
        ].astype(str)

        if n_grids is not None:
            if n_grids < MIN_TRAIN_SITES:
                raise ValueError(
                    f"{source}: n_grids={n_grids} < min {MIN_TRAIN_SITES}"
                )
            if n_grids > len(grids_available):
                raise ValueError(
                    f"{source}: n_grids={n_grids} exceeds available "
                    f"grids ({len(grids_available)})"
                )
            picked_grids = _sample_grids(grids_available, n_grids, seed)
            in_grid = x[GRID_COL].astype(str).isin(picked_grids)
            y, x = y.loc[in_grid], x.loc[in_grid]
            meta["n_grids"] = len(picked_grids)
        else:
            meta["n_grids"] = x[GRID_COL].astype(str).nunique()

        x = x.loc[y.index, bioclim_cols]
        meta["n_points"] = len(y)
    else:
        y = _filter_species_by_occurrence(y, min_occ_grid, "grids (full PA)")

        if n_grids is not None and n_grids < len(y):
            if n_grids < MIN_TRAIN_SITES:
                raise ValueError(
                    f"{source}: n_grids={n_grids} < min {MIN_TRAIN_SITES}"
                )
            pick = _sample_grids(y.index, n_grids, seed)
            y, x = y.loc[pick], x.loc[pick]
        elif n_grids is not None and n_grids > len(y):
            raise ValueError(
                f"{source}: n_grids={n_grids} exceeds available grids ({len(y)})"
            )

        x = x.loc[y.index, bioclim_cols]
        meta["n_grids"] = len(y)
        meta["n_points"] = None

    return x, y, meta


def _load_pac_eval(
    min_occ_grid: int = MIN_SPECIES_LIST_PAC_GRIDS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """PAC eval matrix; default min = catalog threshold (3 cells), not train min (10)."""
    y = pd.read_csv(
        "data/processed/species_matrix_PAC.csv", index_col="site_id"
    )
    x = pd.read_csv(
        "data/processed/covariates_PAC.csv", index_col="site_id"
    )
    common = y.index.intersection(x.index)
    y, x = y.loc[common], x.loc[common]
    y = _filter_species_by_occurrence(y, min_occ_grid, "PAC eval grids")
    x = x.loc[y.index]
    return x, y


def run_model(
    train_source: str,
    size_train: int | None = None,
    seed: int = 42,
    grouped: bool = False,
) -> tuple[pd.Series, int]:
    """
    Train on a source and evaluate per-species AUC on PAC (grid cells).

    size_train: number of COD10X10 cells to include in the effort (for PO,
    all GPS points inside those cells are used for training).

    Returns (per-species AUC on PAC, number of training matrix rows).
    """
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    sfx = file_suffix(train_source, grouped)
    ckpt_tag = (
        f"{train_source}{sfx}_n{size_train if size_train is not None else 'all'}"
    )
    ckpt_path = os.path.join(output_dir, f"best_model_{ckpt_tag}.pt")

    x_train, y_train, meta = _load_source_subset(
        train_source, size_train, seed=seed, grouped=grouped
    )
    x_pac, y_pac = _load_pac_eval()

    mode = "grouped" if grouped else meta["unit"]
    if meta["unit"] == "point":
        print(
            f"{train_source} ({mode}) — {meta['n_grids']} grids, "
            f"{meta['n_points']:,} points, {y_train.shape[1]} species"
        )
    else:
        print(
            f"{train_source} ({mode}) — {meta['n_grids']} grids, "
            f"{y_train.shape[1]} species"
        )
    print(f"PAC — {y_pac.shape[0]} grids, {y_pac.shape[1]} species")

    X_train_bio = x_train.values.astype(np.float32)
    Y_train = y_train.values.astype(np.float32)
    X_pac_bio = x_pac.values.astype(np.float32)
    Y_pac = y_pac.values.astype(np.float32)

    scaler = StandardScaler().fit(X_train_bio)
    X_train = scaler.transform(X_train_bio).astype(np.float32)
    X_pac = scaler.transform(X_pac_bio).astype(np.float32)

    X_train_t = torch.tensor(X_train)
    Y_train_t = torch.tensor(Y_train)

    train_dataset = TensorDataset(X_train_t, Y_train_t)
    n = len(train_dataset)
    if n < 2:
        train_ds = val_ds = train_dataset
    else:
        n_train = max(1, int(0.8 * n))
        n_val = n - n_train
        if n_val == 0:
            n_train, n_val = n - 1, 1
        train_ds, val_ds = random_split(
            train_dataset,
            [n_train, n_val],
            generator=torch.Generator().manual_seed(42),
        )

    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=64)

    input_size = X_train.shape[1]
    train_output_size = Y_train.shape[1]

    model = MLP(
        input_size=input_size,
        hidden_size=200,
        output_size=train_output_size,
        hidden_layers=2,
    )

    if train_source == "PAU":
        criterion = BalancedBCELoss()
    elif train_source in ("POU", "POV"):
        criterion = DeepMaxEntLoss()
    else:
        raise ValueError(f"Invalid train source: {train_source}")

    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience=5, factor=0.5
    )

    def run_epoch(loader, train=True):
        model.train(train)
        total_loss, total = 0.0, 0
        with torch.set_grad_enabled(train):
            for xb, yb in loader:
                logits = model(xb)
                loss = criterion(logits, yb)
                if train:
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                total_loss += loss.item() * len(xb)
                total += len(xb)
        return total_loss / total

    num_epochs = 10
    best_val_loss = float("inf")
    for epoch in range(1, num_epochs + 1):
        train_loss = run_epoch(train_loader, train=True)
        val_loss = run_epoch(val_loader, train=False)
        scheduler.step(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), ckpt_path)

        if epoch % 10 == 0:
            print(
                f"Epoch {epoch:3d} | train loss {train_loss:.4f} "
                f"val loss {val_loss:.4f}"
            )

    print(f"\nBest val loss ({train_source}): {best_val_loss:.4f}  →  {ckpt_path}")

    if not os.path.isfile(ckpt_path):
        raise RuntimeError(
            f"No checkpoint saved. Check NaN in covariates for {train_source}."
        )

    model.load_state_dict(torch.load(ckpt_path, weights_only=True))
    model.eval()

    train_species = list(y_train.columns)
    with torch.no_grad():
        y_score = model(torch.tensor(X_pac)).sigmoid().numpy()

    y_score_df = pd.DataFrame(
        y_score, columns=train_species, index=y_pac.index
    )

    def per_species_auc(
        y_true: pd.DataFrame, y_pred: pd.DataFrame, species: List[str]
    ) -> Dict[str, float]:
        scores: Dict[str, float] = {}
        for sp in species:
            if sp not in y_pred.columns or sp not in y_true.columns:
                scores[sp] = np.nan
                continue
            try:
                auc = roc_auc_score(y_true[sp].values, y_pred[sp].values)
            except ValueError:
                auc = np.nan
            scores[sp] = auc
        return scores

    auc_scores = per_species_auc(y_pac, y_score_df, list(y_pac.columns))
    auc_series = pd.Series(auc_scores).sort_values(ascending=False)
    auc_valid = auc_series.dropna()

    print(
        f"\nPer-species AUC on PAC  (trained on {len(train_species)} species, "
        f"valid AUC n={len(auc_valid)})"
    )
    print(auc_valid.head(20).to_string())
    print(f"\nMean AUC  : {auc_valid.mean():.4f}")
    print(f"Median AUC: {auc_valid.median():.4f}")
    print(f"Skipped   : {auc_series.isna().sum()} species")

    n_tag = meta["n_points"] if meta["n_points"] else meta["n_grids"]
    out_path = os.path.join(
        output_dir,
        f"auc_per_species_PAC_train_{train_source}_n{n_tag}.csv",
    )
    auc_series.to_csv(out_path, header=["AUC"])

    return auc_series, len(y_train)


if __name__ == "__main__":
    auc, n_rows = run_model("PAU", size_train=None)
    print(f"Mean AUC: {auc.mean():.4f}  (n={auc.notna().sum()} species, rows={n_rows})")
