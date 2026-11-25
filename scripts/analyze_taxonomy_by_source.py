#!/usr/bin/env python3
"""Taxonomic composition (family/order) across PAC/PAU/POV/POU.

Uses GBIF backbone family/order from taxonkey_list.xlsx joined to curated keep
list (en_clave_o_aloctonas=SI). Species-level: a taxon counts in a class if
n_obs in that class > 0 (curated_taxon_keep.csv).

Outputs under output/curated/:
  taxonomy_species_by_family_source.csv
  taxonomy_species_by_order_source.csv
  taxonomy_family_enrichment_vs_pac.csv
  taxonomy_by_source_family.png
  taxonomy_by_source_order.png
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plot_style import apply as apply_plot_style  # noqa: E402

CLASSES = ("PAC", "PAU", "POV", "POU")
KEEP_CSV = ROOT / "data" / "raw" / "curated_taxon_keep.csv"
TAXON_XLSX = ROOT / "data" / "raw" / "taxonkey_list.xlsx"
OUT_DIR = ROOT / "output" / "curated"


def _binomial(name: str) -> str:
    s = re.sub(r"^×\s*", "", str(name).strip())
    parts = s.split()
    if len(parts) < 2:
        return s
    if parts[1] in {"×", "x"} and len(parts) >= 3:
        return f"{parts[0]} {parts[2]}"
    return f"{parts[0]} {parts[1]}"


def load_species_taxonomy() -> pd.DataFrame:
    keep = pd.read_csv(KEEP_CSV)
    keep = keep[keep["en_clave_o_aloctonas"] == "SI"].copy()
    tk = pd.read_excel(TAXON_XLSX)
    tk = tk.dropna(subset=["accepted_species_name"]).copy()
    tk["bin"] = tk["accepted_species_name"].map(_binomial)
    tax = (
        tk.groupby("bin", as_index=False)
        .agg(
            accepted_family=("accepted_family", "first"),
            accepted_order=("accepted_order", "first"),
        )
    )
    keep["bin"] = keep["accepted_species_name"].map(_binomial)
    m = keep.merge(tax, on="bin", how="left")
    miss = m["accepted_family"].isna().sum()
    if miss:
        print(f"WARNING: {miss} curated taxa without family in taxonkey_list")
    return m


def species_composition(df: pd.DataFrame, rank: str) -> pd.DataFrame:
    """Share of species per taxonomic rank within each observation class."""
    rows = []
    for cls in CLASSES:
        sub = df[df[cls] > 0]
        n = len(sub)
        vc = sub[rank].fillna("Unknown").value_counts()
        for taxon, count in vc.items():
            rows.append(
                {
                    "class": cls,
                    "rank": rank,
                    "taxon": taxon,
                    "n_species": int(count),
                    "n_species_in_class": n,
                    "share": float(count) / n if n else 0.0,
                }
            )
    return pd.DataFrame(rows)


def enrichment_vs_pac(comp: pd.DataFrame, rank: str, top_n: int = 25) -> pd.DataFrame:
    """log2(share_class / share_PAC) for taxa among top-N by PAC species count."""
    pac = comp[(comp["class"] == "PAC") & (comp["rank"] == rank)].set_index("taxon")
    top = pac.sort_values("n_species", ascending=False).head(top_n).index.tolist()
    rows = []
    eps = 1e-6
    for cls in CLASSES:
        sub = comp[(comp["class"] == cls) & (comp["rank"] == rank)].set_index("taxon")
        for taxon in top:
            p_pac = float(pac.loc[taxon, "share"]) if taxon in pac.index else 0.0
            p = float(sub.loc[taxon, "share"]) if taxon in sub.index else 0.0
            rows.append(
                {
                    "class": cls,
                    "rank": rank,
                    "taxon": taxon,
                    "share": p,
                    "share_pac": p_pac,
                    "log2_enrichment_vs_pac": float(np.log2((p + eps) / (p_pac + eps))),
                    "n_species": int(sub.loc[taxon, "n_species"])
                    if taxon in sub.index
                    else 0,
                }
            )
    return pd.DataFrame(rows)


def _plot_heatmap(
    enrich: pd.DataFrame,
    rank: str,
    out_path: Path,
    title: str,
) -> None:
    apply_plot_style()
    # pivot: taxon × class
    mat = enrich.pivot(index="taxon", columns="class", values="log2_enrichment_vs_pac")
    mat = mat.reindex(columns=list(CLASSES))
    # order taxa by PAC share descending (already top list order in enrich)
    pac_order = (
        enrich[enrich["class"] == "PAC"]
        .sort_values("share", ascending=False)["taxon"]
        .tolist()
    )
    mat = mat.reindex(pac_order)

    fig, ax = plt.subplots(figsize=(5.2, max(4.0, 0.28 * len(mat))))
    vmax = float(np.nanmax(np.abs(mat.values))) if len(mat) else 1.0
    vmax = max(vmax, 0.5)
    im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(CLASSES)))
    ax.set_xticklabels(CLASSES)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index.tolist(), fontsize=7)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r"$\log_2(\mathrm{share}/\mathrm{share}_{PAC})$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"Wrote {out_path}")


def _plot_share_bars(
    comp: pd.DataFrame,
    rank: str,
    out_path: Path,
    top_n: int,
    title: str,
) -> None:
    apply_plot_style()
    pac = comp[(comp["class"] == "PAC") & (comp["rank"] == rank)]
    top = pac.sort_values("n_species", ascending=False).head(top_n)["taxon"].tolist()
    sub = comp[(comp["rank"] == rank) & (comp["taxon"].isin(top))].copy()
    # stacked shares of top taxa + Other
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    bottoms = np.zeros(len(CLASSES))
    # color cycle
    cmap = plt.get_cmap("tab20")
    for i, taxon in enumerate(top):
        shares = []
        for cls in CLASSES:
            row = sub[(sub["class"] == cls) & (sub["taxon"] == taxon)]
            shares.append(float(row["share"].iloc[0]) if len(row) else 0.0)
        ax.bar(
            CLASSES,
            shares,
            bottom=bottoms,
            label=taxon,
            color=cmap(i % 20),
            width=0.7,
        )
        bottoms = bottoms + np.array(shares)
    other = 1.0 - bottoms
    ax.bar(CLASSES, other, bottom=bottoms, label="Other", color="#bbbbbb", width=0.7)
    ax.set_ylabel("Share of species")
    ax.set_ylim(0, 1)
    ax.set_title(title)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=6,
        frameon=False,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = load_species_taxonomy()
    print(
        f"Curated spp={len(df)}  families={df['accepted_family'].nunique()}  "
        f"orders={df['accepted_order'].nunique()}"
    )
    for cls in CLASSES:
        print(f"  {cls}: {(df[cls] > 0).sum()} spp with ≥1 obs")

    fam = species_composition(df, "accepted_family")
    ord_ = species_composition(df, "accepted_order")
    fam.to_csv(args.out_dir / "taxonomy_species_by_family_source.csv", index=False)
    ord_.to_csv(args.out_dir / "taxonomy_species_by_order_source.csv", index=False)

    enr_f = enrichment_vs_pac(fam, "accepted_family", top_n=args.top_n)
    enr_o = enrichment_vs_pac(ord_, "accepted_order", top_n=min(args.top_n, 15))
    enr_f.to_csv(args.out_dir / "taxonomy_family_enrichment_vs_pac.csv", index=False)
    enr_o.to_csv(args.out_dir / "taxonomy_order_enrichment_vs_pac.csv", index=False)

    _plot_heatmap(
        enr_f,
        "accepted_family",
        args.out_dir / "taxonomy_by_source_family.png",
        "Family enrichment vs PAC (species share)",
    )
    _plot_heatmap(
        enr_o,
        "accepted_order",
        args.out_dir / "taxonomy_by_source_order.png",
        "Order enrichment vs PAC (species share)",
    )
    _plot_share_bars(
        fam,
        "accepted_family",
        args.out_dir / "taxonomy_family_share_stacked.png",
        top_n=12,
        title="Species composition by family",
    )
    _plot_share_bars(
        ord_,
        "accepted_order",
        args.out_dir / "taxonomy_order_share_stacked.png",
        top_n=10,
        title="Species composition by order",
    )

    # Brief console: largest |enrichment| for POU/POV vs PAC
    for cls in ("PAU", "POV", "POU"):
        e = enr_f[enr_f["class"] == cls].copy()
        e["abs"] = e["log2_enrichment_vs_pac"].abs()
        top = e.sort_values("abs", ascending=False).head(5)
        print(f"\nTop |enrichment| families in {cls} vs PAC:")
        for _, r in top.iterrows():
            print(
                f"  {r['taxon']}: log2={r['log2_enrichment_vs_pac']:+.2f} "
                f"share={r['share']:.3f} (PAC {r['share_pac']:.3f})"
            )


if __name__ == "__main__":
    main()
