#!/usr/bin/env python3
"""Illustrative figures for Beyond PO and PA (curated Q1–Q3).

Writes PNGs to output/curated/figures/. Numbers only from curated CSVs.
No author names / titles in artwork (integrity).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
CUR = ROOT / "output" / "curated"
OUT = CUR / "figures"

# Palette: levels, not purple-AI default
C_PAC = "#1b4f72"
C_PAU = "#2874a6"
C_POV = "#148f77"
C_POU = "#b9770e"
C_RARE = "#6c3483"
C_MED = "#2471a3"
C_COM = "#117a65"
SOURCE_COLORS = {"PAC": C_PAC, "PAU": C_PAU, "POV": C_POV, "POU": C_POU}


def _style() -> None:
    try:
        from plot_style import apply

        apply()
    except Exception:
        plt.rcParams.update(
            {
                "font.size": 10,
                "axes.spines.top": False,
                "axes.spines.right": False,
                "figure.facecolor": "white",
                "axes.facecolor": "white",
            }
        )


def fig_levels_composition() -> None:
    """Q1: stacked % species by PAC rarity bin across levels."""
    df = pd.read_csv(CUR / "rarity_by_source.csv")
    order = ["PAC", "PAU", "POV", "POU"]
    bins = ["rare", "medium", "common"]
    colors = [C_RARE, C_MED, C_COM]

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    x = np.arange(len(order))
    bottom = np.zeros(len(order))
    for b, col in zip(bins, colors):
        vals = []
        for src in order:
            row = df[
                (df["source"] == src)
                & (df["pac_bin"] == b)
                & (df["mode"].isin(["grid", "point"]))
            ]
            # prefer grid for PAC/PAU, point for PO
            if src in ("PAC", "PAU"):
                row = df[(df["source"] == src) & (df["pac_bin"] == b) & (df["mode"] == "grid")]
            else:
                row = df[(df["source"] == src) & (df["pac_bin"] == b) & (df["mode"] == "point")]
            vals.append(float(row["pct_species"].iloc[0]) if len(row) else 0.0)
        ax.bar(x, vals, bottom=bottom, color=col, width=0.62, label=b, edgecolor="white", lw=0.6)
        bottom += np.array(vals)

    ax.set_xticks(x)
    ax.set_xticklabels(order)
    ax.set_ylabel("% of species (PAC prevalence bins)")
    ax.set_ylim(0, 100)
    ax.legend(title="PAC bin", frameon=False, loc="upper right")
    ax.set_title("Observation levels, not just PO vs PA")
    fig.tight_layout()
    fig.savefig(OUT / "01_levels_composition.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_rarification() -> None:
    """Q1b: PAC-common spp look rarer in which source? (cell grain)."""
    df = pd.read_csv(CUR / "common_pac_in_po_rarification.csv")
    sources = ["PAU", "POV", "POU"]
    data = []
    medians = []
    for s in sources:
        v = df[f"{s}_rarify_ratio"].replace([np.inf, -np.inf], np.nan).dropna()
        v = v.clip(upper=100)  # keep POU tail visible on log scale
        data.append(v.values)
        medians.append(float(np.median(df[f"{s}_rarify_ratio"].replace([np.inf, -np.inf], np.nan).dropna())))

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    parts = ax.violinplot(
        data, positions=range(1, 4), showmeans=False, showmedians=False, showextrema=False
    )
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor([C_PAU, C_POV, C_POU][i])
        body.set_alpha(0.55)
        body.set_edgecolor("black")
        body.set_linewidth(0.6)
    bp = ax.boxplot(data, positions=range(1, 4), widths=0.18, patch_artist=True, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor("white")
        patch.set_alpha(0.9)
    ax.axhline(1.0, color="#7f8c8d", ls="--", lw=1)
    ax.set_yscale("log")
    ax.set_xticks(range(1, 4))
    ax.set_xticklabels(sources)
    ax.set_ylabel("PAC cell % / source cell % (log; clipped ≤100)")
    ax.set_title("Do sources make PAC-common species look rarer?")
    for i, m in enumerate(medians):
        ax.text(i + 1, m * 1.35, f"median {m:.1f}", ha="center", fontsize=8, color="#2c3e50")
    ax.annotate(
        "POU rarifies commons",
        xy=(3, medians[2]),
        xytext=(2.0, 70),
        fontsize=9,
        arrowprops=dict(arrowstyle="->", color="#922b21"),
        color="#922b21",
    )
    fig.tight_layout()
    fig.savefig(OUT / "02_rarification_commons.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_pa_investment() -> None:
    """Q2: AUC gain from adding PAU at fixed POV/POU."""
    mix = pd.read_csv(CUR / "optimal_model_auc_grouped.csv")
    mix = mix[mix["model_type"] == "mix"].copy()

    panels = [
        ("POV100 + POU100", 1.0, 1.0),
        ("POV100 only", 1.0, 0.0),
        ("POU100 only", 0.0, 1.0),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.6), sharey=True)
    for ax, (title, pov, pou) in zip(axes, panels):
        sub = mix[(mix["pct_POV"] == pov) & (mix["pct_POU"] == pou)].sort_values("pct_PAU")
        ax.plot(
            sub["pct_PAU"] * 100,
            sub["mean_auc_pac"],
            "-o",
            color=C_PAU,
            lw=2,
            ms=7,
            markerfacecolor="white",
            markeredgewidth=1.5,
        )
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("% PAU kept")
        ax.set_xticks([0, 50, 100])
        ax.set_xlim(-5, 105)
        ax.grid(True, axis="y", alpha=0.3)
        for _, r in sub.iterrows():
            ax.text(
                r["pct_PAU"] * 100,
                r["mean_auc_pac"] + 0.004,
                f"{r['mean_auc_pac']:.3f}",
                ha="center",
                fontsize=7.5,
            )
    axes[0].set_ylabel("Mean AUC on PAC")
    fig.suptitle("Value of incomplete PA (PAU) at cell grain", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "03_pa_investment.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_inventory_vs_sdm() -> None:
    """Q3: same mixes — species found vs AUC."""
    inv = pd.read_csv(CUR / "species_accumulation_mix.csv")
    sdm = pd.read_csv(CUR / "optimal_model_auc_grouped.csv")
    sdm = sdm[sdm["model_type"] == "mix"][["label", "mean_auc_pac", "pct_PAU", "pct_POV", "pct_POU"]]
    m = inv.merge(sdm, on=["label", "pct_PAU", "pct_POV", "pct_POU"])

    # color by POV level
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    for pov, marker, label in [
        (0.0, "o", "POV 0%"),
        (0.5, "s", "POV 50%"),
        (1.0, "D", "POV 100%"),
    ]:
        sub = m[m["pct_POV"] == pov]
        ax.scatter(
            sub["n_species"],
            sub["mean_auc_pac"],
            s=55 + 40 * sub["pct_PAU"],
            c=[C_POU if pov == 0 else (C_MED if pov == 0.5 else C_POV)] * len(sub),
            marker=marker,
            edgecolors="black",
            linewidths=0.5,
            alpha=0.85,
            label=label,
        )
    # mark optima
    best_inv = m.loc[m["n_species"].idxmax()]
    best_sdm = m.loc[m["mean_auc_pac"].idxmax()]
    ax.scatter(
        [best_inv["n_species"]],
        [best_inv["mean_auc_pac"]],
        s=180,
        facecolors="none",
        edgecolors="#922b21",
        linewidths=2,
        zorder=5,
    )
    ax.annotate(
        "max inventory",
        (best_inv["n_species"], best_inv["mean_auc_pac"]),
        textcoords="offset points",
        xytext=(-50, 12),
        fontsize=8,
        color="#922b21",
    )
    if best_sdm["label"] != best_inv["label"]:
        ax.annotate(
            "max SDM",
            (best_sdm["n_species"], best_sdm["mean_auc_pac"]),
            textcoords="offset points",
            xytext=(8, -14),
            fontsize=8,
            color="#1a5276",
        )
    ax.set_xlabel("Unique species observed (inventory mix)")
    ax.set_ylabel("Mean AUC on PAC (SDM mix)")
    ax.set_title("Optimal mix: inventory vs prediction")
    ax.legend(frameon=False, title="POV share")
    ax.text(
        0.02,
        0.02,
        "Marker size ∝ PAU share",
        transform=ax.transAxes,
        fontsize=8,
        color="#566573",
    )
    fig.tight_layout()
    fig.savefig(OUT / "04_inventory_vs_sdm.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_effort_rare() -> None:
    """Q1/Q2: learning curves for rare species only (clean)."""
    r = pd.read_csv(CUR / "sampling_effort_auc.csv")
    r = r[(r["prevalence_bin"] == "rare") & (r["data_mode"] == "point")]
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for src, c in [("POV", C_POV), ("PAU", C_PAU), ("POU", C_POU)]:
        sub = r[r["train_source"] == src].sort_values("n_grids")
        ax.plot(sub["n_grids"], sub["mean_auc_pac"], "-o", color=c, label=src, lw=2, ms=5)
    ax.set_xscale("log")
    ax.set_xlabel("Training cells sampled (log)")
    ax.set_ylabel("Mean AUC on PAC (rare spp)")
    ax.set_title("Who informs rare species?")
    ax.legend(frameon=False)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "05_effort_rare_point.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_exclusives() -> None:
    """Complementarity PAU vs PO for finding species."""
    ex = pd.read_csv(CUR / "species_source_exclusives.csv", index_col=0)["n"]
    only_pau = int(ex["only_PAU"])
    only_po = int(ex["only_PO"])
    both = int(ex["PAU_and_PO"])

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    # simple three bars
    labs = ["Only PAU", "PAU ∩ PO", "Only PO\n(POV∪POU)"]
    vals = [only_pau, both, only_po]
    cols = [C_PAU, "#5d6d7e", C_POV]
    bars = ax.bar(labs, vals, color=cols, width=0.55, edgecolor="white")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 15, str(v), ha="center", fontsize=11, fontweight="bold")
    ax.set_ylabel("Species")
    ax.set_title("Unique species by observation family")
    ax.set_ylim(0, max(vals) * 1.15)
    fig.tight_layout()
    fig.savefig(OUT / "06_species_exclusives.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_levels_schematic() -> None:
    """Conceptual strip: four levels on one axis of quality/structure."""
    fig, ax = plt.subplots(figsize=(8.2, 2.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis("off")

    levels = [
        (1.2, "POU", "PO unvalidated", C_POU, "GBIF points"),
        (3.5, "POV", "PO validated", C_POV, "Biodibal"),
        (5.8, "PAU", "PA incomplete", C_PAU, "partial atlas"),
        (8.1, "PAC", "PA complete", C_PAC, "eval anchor"),
    ]
    ax.annotate(
        "",
        xy=(9.3, 1.5),
        xytext=(0.5, 1.5),
        arrowprops=dict(arrowstyle="->", color="#34495e", lw=1.4),
    )
    ax.text(5, 0.35, "structure / validation of observation", ha="center", fontsize=9, color="#566573")
    for x, code, name, col, note in levels:
        box = FancyBboxPatch(
            (x - 0.85, 1.05),
            1.7,
            1.15,
            boxstyle="round,pad=0.02,rounding_size=0.15",
            facecolor=col,
            edgecolor="white",
            alpha=0.9,
            linewidth=1.5,
        )
        ax.add_patch(box)
        ax.text(x, 1.75, code, ha="center", va="center", color="white", fontsize=12, fontweight="bold")
        ax.text(x, 1.35, name, ha="center", va="center", color="white", fontsize=7.5)
        ax.text(x, 0.85, note, ha="center", va="top", fontsize=7.5, color="#2c3e50")
    ax.set_title("Beyond PO and PA: four observation-quality levels", fontsize=11, pad=8)
    fig.tight_layout()
    fig.savefig(OUT / "00_levels_schematic.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig_top_mixes() -> None:
    """Horizontal top-10 SDM mixes (readable 2D; not a fake 4D cube)."""
    r = pd.read_csv(CUR / "optimal_model_auc_grouped.csv")
    r = r[r["model_type"] == "mix"].sort_values("mean_auc_pac", ascending=True).tail(10)
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    colors = [C_POV if p == 1 else (C_MED if p == 0.5 else C_POU) for p in r["pct_POV"]]
    ax.barh(r["label"], r["mean_auc_pac"], color=colors, edgecolor="white")
    ax.set_xlabel("Mean AUC on PAC")
    ax.set_xlim(0.68, 0.80)
    ax.set_title("Top SDM mixes (cell grain, curated)")
    for y, v in enumerate(r["mean_auc_pac"]):
        ax.text(v + 0.002, y, f"{v:.3f}", va="center", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "07_top_sdm_mixes.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    alt = OUT / "07_top_sdm_mixes_alt_view.png"
    if alt.exists():
        alt.unlink()


def main() -> int:
    _style()
    OUT.mkdir(parents=True, exist_ok=True)
    fig_levels_schematic()
    fig_levels_composition()
    fig_rarification()
    fig_pa_investment()
    fig_inventory_vs_sdm()
    fig_effort_rare()
    fig_exclusives()
    fig_top_mixes()
    print(f"Wrote figures → {OUT}")
    for p in sorted(OUT.glob("*.png")):
        print(f"  {p.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
