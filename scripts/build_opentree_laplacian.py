#!/usr/bin/env python3
"""Build phylogenetic Laplacian from OpenTree for curated species columns.

Caches:
  data/processed/curated/opentree/match_report.csv
  data/processed/curated/opentree/induced_subtree.newick
  data/processed/curated/opentree/laplacian.npz  (L, species, sigma, n_matched)
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from preprocessing import (  # noqa: E402
    set_processed_dir,
    species_matrix_path,
)

OT_TNRS = "https://api.opentreeoflife.org/v3/tnrs/match_names"
OT_INDUCED = "https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree"
CONTEXT = "Land plants"
TNRS_CHUNK = 100
REQUEST_PAUSE_S = 0.25


def _binomial(name: str) -> str:
    """Genus + epithet; drop hybrid marks / authorities."""
    s = re.sub(r"^×\s*", "", name.strip())
    parts = s.split()
    if len(parts) < 2:
        return s
    # skip hybrid marker as second token
    if parts[1] in {"×", "x"} and len(parts) >= 3:
        return f"{parts[0]} {parts[2]}"
    return f"{parts[0]} {parts[1]}"


def _tnrs_batch(names: list[str], session: requests.Session) -> list[dict]:
    r = session.post(
        OT_TNRS,
        json={"names": names, "context_name": CONTEXT, "do_approximate_matching": True},
        timeout=120,
    )
    r.raise_for_status()
    return r.json().get("results", [])


def match_species(species: list[str]) -> pd.DataFrame:
    """Match catalogue names to OTT ids. Prefer exact, else binomial query."""
    session = requests.Session()
    rows: list[dict] = []

    # Query both full string and binomial (dedupe queries).
    query_for: dict[str, list[str]] = {}
    for sp in species:
        for q in (sp, _binomial(sp)):
            query_for.setdefault(q, []).append(sp)

    queries = list(query_for.keys())
    ott_by_query: dict[str, tuple[int | None, str, str]] = {}

    for i in range(0, len(queries), TNRS_CHUNK):
        chunk = queries[i : i + TNRS_CHUNK]
        results = _tnrs_batch(chunk, session)
        time.sleep(REQUEST_PAUSE_S)
        by_name = {res["name"]: res for res in results}
        for q in chunk:
            res = by_name.get(q)
            if not res or not res.get("matches"):
                ott_by_query[q] = (None, "no_match", "")
                continue
            # Prefer taxonomic match score / exact
            matches = res["matches"]
            matches_sorted = sorted(
                matches,
                key=lambda m: (
                    0 if m.get("is_synonym") else 1,
                    float(m.get("score", 0.0)),
                ),
                reverse=True,
            )
            best = matches_sorted[0]
            taxon = best.get("taxon") or {}
            ott = taxon.get("ott_id")
            unique = taxon.get("unique_name") or taxon.get("name") or ""
            flag = "matched"
            if best.get("is_synonym"):
                flag = "synonym"
            ott_by_query[q] = (int(ott) if ott is not None else None, flag, unique)

    for sp in species:
        full = ott_by_query.get(sp, (None, "no_match", ""))
        binom = ott_by_query.get(_binomial(sp), (None, "no_match", ""))
        # Prefer full-name match when both exist
        if full[0] is not None:
            ott, flag, unique = full
            query_used = sp
        elif binom[0] is not None:
            ott, flag, unique = binom
            query_used = _binomial(sp)
            flag = f"binomial_{flag}"
        else:
            ott, flag, unique, query_used = None, "no_match", "", sp
        rows.append(
            {
                "species": sp,
                "query": query_used,
                "ott_id": ott,
                "match_flag": flag,
                "ott_unique_name": unique,
            }
        )

    return pd.DataFrame(rows)


def fetch_induced_subtree(ott_ids: list[int]) -> str:
    """Induced subtree; drop pruned/unknown OTT ids and retry."""
    session = requests.Session()
    remaining = list(ott_ids)
    last_detail: object = None
    for attempt in range(8):
        r = session.post(
            OT_INDUCED,
            json={"ott_ids": remaining, "label_format": "name_and_id"},
            timeout=300,
        )
        if r.status_code < 400:
            data = r.json()
            newick = data.get("newick") or data.get("subtree")
            if not newick:
                raise RuntimeError(
                    f"No newick in induced_subtree response keys={list(data)}"
                )
            if attempt:
                print(
                    f"  induced_subtree OK after dropping pruned "
                    f"({len(ott_ids) - len(remaining)} ids removed)"
                )
            return newick
        try:
            detail = r.json()
        except Exception:
            detail = r.text
        last_detail = detail
        unknown = {}
        if isinstance(detail, dict):
            unknown = detail.get("unknown") or {}
        drop = set()
        for key in unknown:
            m = re.search(r"(\d+)", str(key))
            if m:
                drop.add(int(m.group(1)))
        # Also parse message "node_id 'ottN' was not found"
        if isinstance(detail, dict):
            msg = str(detail.get("message", ""))
            for m in re.finditer(r"ott(\d+)", msg):
                drop.add(int(m.group(1)))
        if not drop:
            raise RuntimeError(
                f"induced_subtree failed ({r.status_code}): {detail}"
            )
        before = len(remaining)
        remaining = [o for o in remaining if o not in drop]
        print(
            f"  induced_subtree attempt {attempt + 1}: drop {before - len(remaining)} "
            f"unknown/pruned OTT; left {len(remaining)}"
        )
        if len(remaining) < 2:
            break
        time.sleep(REQUEST_PAUSE_S)
    raise RuntimeError(f"induced_subtree failed after retries: {last_detail}")


def _parse_ott_from_label(label: str) -> int | None:
    m = re.search(r"ott(\d+)", label.replace(" ", "_"))
    if m:
        return int(m.group(1))
    m = re.search(r"_ott_?(\d+)$", label)
    if m:
        return int(m.group(1))
    return None


def cophenetic_from_newick(newick: str) -> tuple[np.ndarray, list[int]]:
    import dendropy

    tree = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    # OpenTree induced trees often lack lengths; default unit edges.
    for edge in tree.preorder_edge_iter():
        if edge.length is None:
            edge.length = 1.0

    taxa = []
    otts: list[int] = []
    for leaf in tree.leaf_node_iter():
        if leaf.taxon is None:
            continue
        label = str(leaf.taxon.label)
        ott = _parse_ott_from_label(label)
        if ott is None:
            continue
        taxa.append(leaf.taxon)
        otts.append(ott)

    n = len(taxa)
    D = np.zeros((n, n), dtype=np.float64)
    pdm = tree.phylogenetic_distance_matrix()
    for i in range(n):
        for j in range(i + 1, n):
            d = pdm.patristic_distance(taxa[i], taxa[j])
            D[i, j] = D[j, i] = float(d)
    return D, otts


def build_laplacian(
    species: list[str],
    match_df: pd.DataFrame,
    D_ott: np.ndarray,
    ott_order: list[int],
) -> tuple[np.ndarray, float, int]:
    """Align L to species list; unmatched / missing tips isolated."""
    k = len(species)
    L = np.zeros((k, k), dtype=np.float64)
    ott_to_i = {ott: i for i, ott in enumerate(ott_order)}

    matched_idx: list[int] = []
    matched_ott_pos: list[int] = []
    for si, sp in enumerate(species):
        row = match_df.loc[match_df["species"] == sp].iloc[0]
        ott = row["ott_id"]
        if pd.isna(ott):
            continue
        ott = int(ott)
        if ott not in ott_to_i:
            continue
        matched_idx.append(si)
        matched_ott_pos.append(ott_to_i[ott])

    n_matched = len(matched_idx)
    if n_matched < 2:
        return L, float("nan"), n_matched

    # Submatrix of distances for matched tips (species order of matched_idx)
    sub = D_ott[np.ix_(matched_ott_pos, matched_ott_pos)]
    # Collapse duplicate OTT (synonyms mapping to same tip): keep first
    # already unique ott_pos; if two species share ott, they get identical rows
    iu = np.triu_indices(n_matched, k=1)
    dists = sub[iu]
    dists = dists[np.isfinite(dists) & (dists > 0)]
    sigma = float(np.median(dists)) if len(dists) else 1.0
    if sigma <= 0:
        sigma = 1.0

    A_sub = np.exp(-sub / sigma)
    np.fill_diagonal(A_sub, 0.0)

    A = np.zeros((k, k), dtype=np.float64)
    for a, ia in enumerate(matched_idx):
        for b, ib in enumerate(matched_idx):
            A[ia, ib] = A_sub[a, b]

    deg = A.sum(axis=1)
    L = np.diag(deg) - A
    return L, sigma, n_matched


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--processed-dir",
        default="data/processed/curated",
        help="Matrices dir (default curated)",
    )
    parser.add_argument(
        "--source",
        default="PAU",
        help="Species columns from this matrix (default PAU)",
    )
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Cache dir (default <processed-dir>/opentree)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rebuild even if cache exists",
    )
    args = parser.parse_args()

    set_processed_dir(args.processed_dir)
    out_dir = Path(args.out_dir or os.path.join(args.processed_dir, "opentree"))
    out_dir.mkdir(parents=True, exist_ok=True)

    lap_path = out_dir / "laplacian.npz"
    match_path = out_dir / "match_report.csv"
    nwk_path = out_dir / "induced_subtree.newick"

    y = pd.read_csv(species_matrix_path(args.source, grouped=False), index_col=0, nrows=0)
    species = list(y.columns)
    print(f"Species columns from {args.source}: {len(species)}")

    if match_path.is_file() and not args.force:
        match_df = pd.read_csv(match_path)
        print(f"Reusing match report {match_path}")
    else:
        print("TNRS matching via OpenTree…")
        match_df = match_species(species)
        match_df.to_csv(match_path, index=False)

    n_hit = match_df["ott_id"].notna().sum()
    rate = n_hit / max(len(species), 1)
    print(f"TNRS match: {n_hit}/{len(species)} ({rate:.1%})")

    ott_ids = sorted({int(x) for x in match_df["ott_id"].dropna().unique()})
    print(f"Unique OTT ids: {len(ott_ids)}")

    if nwk_path.is_file() and not args.force:
        newick = nwk_path.read_text().strip()
        print(f"Reusing newick {nwk_path}")
    else:
        print("Fetching induced subtree…")
        newick = fetch_induced_subtree(ott_ids)
        nwk_path.write_text(newick + "\n")

    if lap_path.is_file() and not args.force:
        print(f"Cache exists {lap_path}; use --force to rebuild")
        z = np.load(lap_path, allow_pickle=True)
        print(
            f"Cached L shape={z['L'].shape} sigma={z['sigma']} "
            f"n_matched={z['n_matched']}"
        )
        return

    print("Cophenetic distances…")
    D_ott, ott_order = cophenetic_from_newick(newick)
    print(f"Tips with OTT in tree: {len(ott_order)}")

    L, sigma, n_matched = build_laplacian(species, match_df, D_ott, ott_order)
    np.savez_compressed(
        lap_path,
        L=L.astype(np.float32),
        species=np.array(species, dtype=object),
        sigma=np.array(sigma),
        n_matched=np.array(n_matched),
        match_rate=np.array(rate),
    )
    print(
        f"Wrote {lap_path}  L={L.shape}  sigma={sigma:.4g}  "
        f"n_matched={n_matched}  match_rate={rate:.1%}"
    )


if __name__ == "__main__":
    main()
