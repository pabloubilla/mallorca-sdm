#!/usr/bin/env python3
"""Materialize flora∪alóctonas keep list and curated occurrence products.

Keep rule (product/decisions.md): en_clave_o_aloctonas == SI in
dataset_clave_aloctonas_review.xlsx#revision. Does not apply rescate_texto.

Writes:
  data/raw/curated_taxon_keep.csv
  data/raw/full_data_saez_curated.gpkg
  data/raw/curated/pa_occurrences.csv
  data/raw/curated/extracted_data.csv
  data/raw/curated/filter_log.txt
"""

from __future__ import annotations

import gzip
import shutil
import sys
from pathlib import Path

import geopandas as gpd
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REVIEW = ROOT / "data" / "raw" / "dataset_clave_aloctonas_review.xlsx"
GPKG_IN = ROOT / "data" / "raw" / "full_data_saez.gpkg"
GPKG_OUT = ROOT / "data" / "raw" / "full_data_saez_curated.gpkg"
KEEP_CSV = ROOT / "data" / "raw" / "curated_taxon_keep.csv"
CURATED_RAW = ROOT / "data" / "raw" / "curated"
PA_OCC_IN = ROOT / "data" / "raw" / "pa_occurrences.csv"
EXTRACTED_IN = ROOT / "data" / "raw" / "extracted_data.csv"
LOG_PATH = CURATED_RAW / "filter_log.txt"


def load_keep(review_path: Path) -> pd.DataFrame:
    rev = pd.read_excel(review_path, sheet_name="revision")
    keep = rev[rev["en_clave_o_aloctonas"].astype(str).str.upper() == "SI"].copy()
    keep["accepted_species_key"] = pd.to_numeric(
        keep["accepted_species_key"], errors="coerce"
    )
    return keep


def _class_counts(df: pd.DataFrame, name_col: str = "scientificName") -> pd.DataFrame:
    if "class" not in df.columns:
        return pd.DataFrame()
    rows = []
    for cls, sub in df.groupby("class"):
        rows.append(
            {
                "class": cls,
                "n_obs": len(sub),
                "n_spp": sub[name_col].nunique(),
            }
        )
    return pd.DataFrame(rows).sort_values("class")


def main() -> int:
    if not REVIEW.is_file():
        print(f"Missing {REVIEW}", file=sys.stderr)
        return 1
    if not GPKG_IN.is_file():
        print(f"Missing {GPKG_IN}", file=sys.stderr)
        return 1

    keep = load_keep(REVIEW)
    keep_names = set(keep["accepted_species_name"].dropna().astype(str))
    keep_out = keep[
        [
            "accepted_species_key",
            "accepted_species_name",
            "canonical",
            "en_clave",
            "nombre_en_clave",
            "en_aloctonas",
            "nombre_en_aloctonas",
            "en_clave_o_aloctonas",
            "n_obs_gpkg",
            "PAC",
            "PAU",
            "POV",
            "POU",
        ]
    ].sort_values("accepted_species_name")
    KEEP_CSV.parent.mkdir(parents=True, exist_ok=True)
    keep_out.to_csv(KEEP_CSV, index=False)
    print(f"Wrote {KEEP_CSV} ({len(keep_out)} taxa)")

    print(f"Reading {GPKG_IN} …")
    gdf = gpd.read_file(GPKG_IN)
    before = _class_counts(gdf)
    n_before_obs, n_before_spp = len(gdf), gdf["scientificName"].nunique()
    gdf_c = gdf[gdf["scientificName"].astype(str).isin(keep_names)].copy()
    after = _class_counts(gdf_c)
    n_after_obs, n_after_spp = len(gdf_c), gdf_c["scientificName"].nunique()
    if GPKG_OUT.exists():
        GPKG_OUT.unlink()
    gdf_c.to_file(GPKG_OUT, driver="GPKG")
    print(
        f"Wrote {GPKG_OUT}: {n_before_obs:,}→{n_after_obs:,} obs; "
        f"{n_before_spp}→{n_after_spp} spp"
    )

    CURATED_RAW.mkdir(parents=True, exist_ok=True)
    lines = [
        f"keep_taxa={len(keep_names)}",
        f"gpkg_obs_before={n_before_obs}",
        f"gpkg_obs_after={n_after_obs}",
        f"gpkg_spp_before={n_before_spp}",
        f"gpkg_spp_after={n_after_spp}",
        "",
        "gpkg_by_class_before:",
        before.to_string(index=False) if not before.empty else "(none)",
        "",
        "gpkg_by_class_after:",
        after.to_string(index=False) if not after.empty else "(none)",
        "",
    ]

    if PA_OCC_IN.is_file():
        pa = pd.read_csv(PA_OCC_IN)
        pa_b = len(pa)
        pa_c = pa[pa["scientificName"].astype(str).isin(keep_names)].copy()
        pa_c.to_csv(CURATED_RAW / "pa_occurrences.csv", index=False)
        msg = (
            f"pa_occurrences: {pa_b:,}→{len(pa_c):,} rows; "
            f"spp {pa['scientificName'].nunique()}→{pa_c['scientificName'].nunique()}"
        )
        lines.append(msg)
        print(msg)
    else:
        lines.append(f"pa_occurrences: MISSING {PA_OCC_IN}")

    if EXTRACTED_IN.is_file():
        # cell_pixel rows have no scientificName (labels from pa_occurrences);
        # filter only point rows.
        parts = []
        n_in = n_out = 0
        for chunk in pd.read_csv(
            EXTRACTED_IN, chunksize=80_000, low_memory=False
        ):
            n_in += len(chunk)
            if "extract_mode" in chunk.columns:
                pix = chunk[chunk["extract_mode"] == "cell_pixel"]
                pts = chunk[chunk["extract_mode"] != "cell_pixel"]
                pts = pts[pts["scientificName"].astype(str).isin(keep_names)]
                sub = pd.concat([pix, pts], ignore_index=True)
            else:
                sub = chunk[chunk["scientificName"].astype(str).isin(keep_names)]
            n_out += len(sub)
            parts.append(sub)
        ex_c = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
        out_ex = CURATED_RAW / "extracted_data.csv"
        out_gz = CURATED_RAW / "extracted_data.csv.gz"
        ex_c.to_csv(out_ex, index=False)
        with open(out_ex, "rb") as src, gzip.open(out_gz, "wb") as dst:
            shutil.copyfileobj(src, dst)
        msg = (
            f"extracted_data: {n_in:,}→{n_out:,} rows "
            f"(wrote {out_ex.name} + {out_gz.name})"
        )
        lines.append(msg)
        print(msg)
    else:
        lines.append(f"extracted_data: MISSING {EXTRACTED_IN}")

    LOG_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {LOG_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
