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

def main():

    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)

    #### For now this script just runs one source and tests in the other
    # PAU -> PAC
    # It should be better to have it be the source as an inputs
    # Keep in mind that depending on the source we might want to use different losses
    # For PO Normally DeepMaxent, for PA use BCE or Balanced BCE 

    # ── Load species matrices and covariates ───────────────────────────────────────────
    y_pau  = pd.read_csv("data/processed/species_matrix_PAU.csv", index_col="site_id")
    y_pac  = pd.read_csv("data/processed/species_matrix_PAC.csv", index_col="site_id")
    x_pau = pd.read_csv("data/processed/covariates_PAU.csv", index_col="site_id")
    x_pac = pd.read_csv("data/processed/covariates_PAC.csv", index_col="site_id")

    print(f"PAU — {y_pau.shape[0]} sites, {y_pau.shape[1]} species")
    print(f"PAC — {y_pac.shape[0]} sites, {y_pac.shape[1]} species")

    # ── Build inputs ──────────────────────────────────────────────────────────────
    # PAU: only Bioclim (we don't have PAC presence at train time)
    X_pau_bio = x_pau.values.astype(np.float32)
    Y_pau     = y_pau.values.astype(np.float32)

    # PAC: only Bioclim (same feature space → model is transferable)
    X_pac_bio = x_pac.values.astype(np.float32)
    Y_pac     = y_pac.values.astype(np.float32)

    # Fit scaler on PAU, apply to both (PAC is held-out — no data leakage)
    scaler = StandardScaler().fit(X_pau_bio)
    X_pau  = scaler.transform(X_pau_bio).astype(np.float32)
    X_pac  = scaler.transform(X_pac_bio).astype(np.float32)

    # ── PAU train / val split (internal, for early stopping only) ─────────────────
    X_pau_t = torch.tensor(X_pau)
    Y_pau_t = torch.tensor(Y_pau)

    pau_dataset = TensorDataset(X_pau_t, Y_pau_t)
    n_train     = int(0.8 * len(pau_dataset))
    n_val       = len(pau_dataset) - n_train
    train_ds, val_ds = random_split(
        pau_dataset, [n_train, n_val],
        generator=torch.Generator().manual_seed(42)
    )

    train_loader = DataLoader(train_ds, batch_size=64, shuffle=True)
    val_loader   = DataLoader(val_ds,   batch_size=64)

    # ── Model ─────────────────────────────────────────────────────────────────────
    # Input: Bioclim only — output: PAU species
    # (head will be replaced for PAC evaluation)
    input_size      = X_pau.shape[1]   # n_bioclim
    pau_output_size = Y_pau.shape[1]   # n_pau species
    pac_output_size = Y_pac.shape[1]   # n_pac species

    model = MLP(
        input_size    = input_size,
        hidden_size   = 200,
        output_size   = pau_output_size,
        hidden_layers = 2,
    )

    criterion = BalancedBCELoss() # nn.BCEWithLogitsLoss() # DeepMaxEntLoss() # BalancedBCELoss() # nn.BCEWithLogitsLoss()
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=5, factor=0.5)

    # ── Training on PAU ───────────────────────────────────────────────────────────
    def run_epoch(loader, train=True):
        model.train(train)
        total_loss, correct, total = 0.0, 0, 0
        with torch.set_grad_enabled(train):
            for xb, yb in loader:
                logits = model(xb)
                loss   = criterion(logits, yb)
                if train:
                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()
                total_loss += loss.item() * len(xb)
                total   += len(xb)
        return total_loss / total

    num_epochs = 10
    best_val_loss = float("inf")
    for epoch in range(1, num_epochs + 1):
        train_loss = run_epoch(train_loader, train=True)
        val_loss   = run_epoch(val_loader,   train=False)
        scheduler.step(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), os.path.join(output_dir, "best_model_pau.pt"))

        if epoch % 10 == 0:
            print(f"Epoch {epoch:3d} | "
                f"train loss {train_loss:.4f} "
                f"val loss {val_loss:.4f}")

    print(f"\nBest val loss (PAU): {best_val_loss:.4f}  →  best_model_pau.pt")



    # ── AUC evaluation on PAC test split ─────────────────────────────────────────
    model.load_state_dict(torch.load(os.path.join(output_dir, "best_model_pau.pt")))
    model.eval()

    # pac_test_loader = DataLoader(TensorDataset(torch.tensor(X_pac), torch.tensor(Y_pac)), batch_size=64)

    # all_logits, all_labels = [], []
    # with torch.no_grad():
    #     for xb, yb in pac_test_loader:
    #         all_logits.append(model(xb).sigmoid().numpy())
    #         all_labels.append(yb.numpy())
    with torch.no_grad():
        y_score = model(torch.tensor(X_pac)).sigmoid().numpy()
    y_true = Y_pac

    y_true_df = pd.DataFrame(y_true, columns=y_pac.columns)


    def per_species_auc(y_true: pd.DataFrame, y_score: np.ndarray, species: List[str]) -> Dict[str, float]:
        scores: Dict[str, float] = {}
        for i, sp in enumerate(species):
            try:
                auc = roc_auc_score(y_true[sp].values, y_score[:, i])
            except ValueError:   # single class in test split
                auc = np.nan
            scores[sp] = auc
        return scores


    auc_scores = per_species_auc(y_true_df, y_score, list(y_pac.columns))
    auc_series = pd.Series(auc_scores).dropna().sort_values(ascending=False)

    print(f"\nPer-species AUC on PAC test set  (n={len(auc_series)} species)")
    print(auc_series.to_string())
    print(f"\nMean AUC  : {auc_series.mean():.4f}")
    print(f"Median AUC: {auc_series.median():.4f}")
    print(f"Skipped   : {sum(np.isnan(v) for v in auc_scores.values())} species (single class in test)")

    auc_series.to_csv(os.path.join(output_dir, "auc_per_species_PAC.csv"), header=["AUC"])


if __name__ == "__main__":
    main()