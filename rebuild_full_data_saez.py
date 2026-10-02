"""Rebuild taxonkey_list + full_data_saez.gpkg (no covariates) incl. Biel.

Mirrors data_preparation.R merge + GBIF name unify. Run with network for
api.gbif.org. Uses portfolio/research venv (geopandas, openpyxl).
"""
from __future__ import annotations

import math
import re
import time
from pathlib import Path
from urllib.parse import quote

import geopandas as gpd
import pandas as pd
import requests
from shapely.geometry import Point

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"

UNC_1KM = 1000 * math.sqrt(2) / 2
UNC_500M = 500 * math.sqrt(2) / 2
CRS = "EPSG:25831"
MERGE_COLS = [
    "eventDate",
    "scientificName",
    "geometry",
    "class",
    "coordinateUncertaintyInMeters",
    "source",
]
FLORA_UTM_RE = re.compile(
    r"^31S\s+([A-Z]{2})\s+(\d+)\s+(\d+)(?:\s+([A-Da-d]))?$", re.I
)
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "mallorca-sdm/rebuild_full_data_saez"})


def normalize_scientific_name(x) -> str | None:
    if pd.isna(x):
        return None
    return re.sub(r"\s+", " ", str(x)).strip()


def normalize_event_date(series: pd.Series) -> pd.Series:
    return pd.to_datetime(series, errors="coerce", utc=True)


def utm_quad_to_code(values: pd.Series) -> pd.Series:
    """\"DD 459 4373\" → \"DD5973\" (Quadricules_GBO / Flora_net)."""
    out = pd.Series([None] * len(values), index=values.index, dtype=object)
    for i, raw in values.items():
        if pd.isna(raw):
            continue
        s = str(raw).strip()
        if not s:
            continue
        parts = s.split()
        if len(parts) < 3:
            continue
        e, n = parts[1], parts[2]
        out.at[i] = f"{parts[0]}{e[-2:]}{n[-2:]}"
    return out


def parse_flora_utm_1x1(values: pd.Series) -> pd.DataFrame:
    rows = []
    for raw in values.astype(str).str.strip():
        m = FLORA_UTM_RE.match(raw)
        if not m:
            rows.append((None, None, None, None))
            continue
        letters, e, n, q = m.group(1), m.group(2), m.group(3), m.group(4)
        code = f"{letters}{e[-2:]}{n[-2:]}"
        q = q.upper() if q else None
        rows.append((e, n, q, code))
    return pd.DataFrame(rows, columns=["utm_e", "utm_n", "quad_500", "UTMCODE1X1"])


def offset_500_xy(quad: pd.Series) -> tuple[pd.Series, pd.Series]:
    q = quad.fillna("").astype(str).str.upper()
    dx = q.map({"A": -250, "B": 250, "C": -250, "D": 250}).fillna(0)
    dy = q.map({"A": -250, "B": -250, "C": 250, "D": 250}).fillna(0)
    return dx.astype(float), dy.astype(float)


def ensure_geom(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gdf = gdf.copy()
    if gdf.geometry.name != "geometry":
        gdf = gdf.set_geometry(gdf.geometry.name)
        gdf = gdf.rename_geometry("geometry")
    return gdf


def load_biodibal() -> gpd.GeoDataFrame:
    gdf = gpd.read_file(RAW / "Biodibal_terrestrial_25831.gpkg")
    gdf = ensure_geom(gdf)
    gdf["coordinateUncertaintyInMeters"] = 1.0
    gdf["eventDate"] = normalize_event_date(gdf["eventDate"])
    gdf["source"] = "Biodibal"
    gdf["class"] = "POV"
    return gdf[MERGE_COLS]


def load_saez() -> gpd.GeoDataFrame:
    df = pd.read_csv(RAW / "Puig_Major.csv", sep=";")
    long = df.melt(
        id_vars=[c for c in df.columns if c not in ("31SDE8206", "31SDE8207")],
        value_vars=["31SDE8206", "31SDE8207"],
        var_name="quad",
        value_name="presence",
    )
    long = long[long["presence"] == 1].copy()
    xy = {
        "31SDE8206": (482507, 4406457),
        "31SDE8207": (483512, 4407435),
    }
    long["XCENTROIDE"] = long["quad"].map(lambda q: xy[q][0])
    long["YCENTROIDE"] = long["quad"].map(lambda q: xy[q][1])
    gdf = gpd.GeoDataFrame(
        {
            "eventDate": pd.NaT,
            "scientificName": long["SP"].values,
            "coordinateUncertaintyInMeters": UNC_500M,
            "source": "Saez",
            "class": "PAC",
        },
        geometry=[Point(xy) for xy in zip(long["XCENTROIDE"], long["YCENTROIDE"])],
        crs=CRS,
    )
    return gdf[MERGE_COLS]


def load_gbif_and_flora(utm_net: gpd.GeoDataFrame):
    gbif = gpd.read_file(RAW / "GBIF_Mallorca.gpkg")
    gbif = ensure_geom(gbif)
    gbif = gbif[gbif["taxonRank"].isin(["SPECIES", "SUBSPECIES", "VARIETY"])]
    gbif = gbif[
        gbif["coordinateUncertaintyInMeters"].notna()
        & (gbif["coordinateUncertaintyInMeters"] < 1000)
    ]
    potential = gbif[gbif["coordinateUncertaintyInMeters"] == 707].copy()
    potential["eventDate"] = normalize_event_date(potential["eventDate"])
    potential["source"] = "GBIF"
    potential["coordinateUncertaintyInMeters"] = UNC_1KM
    potential_bind = potential[
        ["eventDate", "scientificName", "geometry", "source", "coordinateUncertaintyInMeters"]
    ]

    gbif_def = gbif[
        (gbif["coordinateUncertaintyInMeters"] != 707)
        & (gbif["collectionCode"] != "FV-MALLORCA")
    ].copy()
    gbif_def["eventDate"] = normalize_event_date(gbif_def["eventDate"])
    gbif_def["source"] = "GBIF"
    gbif_def["class"] = "POU"
    gbif_po = gbif_def[
        ["eventDate", "scientificName", "geometry", "class", "coordinateUncertaintyInMeters", "source"]
    ]

    saez = load_saez()
    saez_bind = saez[
        ["eventDate", "scientificName", "geometry", "source", "coordinateUncertaintyInMeters"]
    ]
    full_flora = pd.concat([potential_bind, saez_bind], ignore_index=True)
    full_flora = gpd.GeoDataFrame(full_flora, geometry="geometry", crs=CRS)

    flora_net = utm_net[utm_net["Flora"] == 1][["Flora", "geometry"]].copy()
    joined = gpd.sjoin(full_flora, flora_net, how="left", predicate="intersects")
    joined = joined[~joined.index.duplicated(keep="first")]
    flora_pac = joined[joined["Flora"] == 1][
        ["eventDate", "scientificName", "geometry", "source", "coordinateUncertaintyInMeters"]
    ].copy()
    flora_pac["class"] = "PAC"

    flora_pau = gpd.sjoin(
        gpd.GeoDataFrame(potential_bind, geometry="geometry", crs=CRS),
        flora_net,
        how="left",
        predicate="intersects",
    )
    flora_pau = flora_pau[~flora_pau.index.duplicated(keep="first")]
    flora_pau = flora_pau[flora_pau["Flora"].isna()][
        ["eventDate", "scientificName", "geometry", "source", "coordinateUncertaintyInMeters"]
    ].copy()
    flora_pau["class"] = "PAU"
    return flora_pac[MERGE_COLS], flora_pau[MERGE_COLS], gbif_po[MERGE_COLS]


def load_biel(utm_net: gpd.GeoDataFrame) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    quad_path = RAW / "Biel" / "20260620_Quadrícules_GBO.xlsx"
    if not quad_path.exists():
        quad_path = RAW / "Quadricules_GBO.xlsx"
    biel_q = pd.read_excel(quad_path)
    if "UTMCODE1X1" not in biel_q.columns:
        if "UTM" not in biel_q.columns:
            raise ValueError("Biel cuadrículas missing UTMCODE1X1 and UTM")
        biel_q = biel_q.copy()
        biel_q["UTMCODE1X1"] = utm_quad_to_code(biel_q["UTM"])
    pac_codes = set(
        biel_q.loc[(biel_q["Illa"] == "Mallorca") & biel_q["UTMCODE1X1"].notna(), "UTMCODE1X1"]
    )

    # metadata gpkg (cuadrículas joined to net)
    biel_utm = utm_net.merge(biel_q, on="UTMCODE1X1", how="left")
    biel_utm.to_file(RAW / "Biel_data_utm.gpkg", driver="GPKG")

    flora_raw = pd.read_excel(RAW / "Biel" / "20260823_FLORA_GBO_Reduida.xlsx")
    flora_raw = flora_raw[
        (flora_raw["Illa"] == "Mallorca") & flora_raw["UTM 1x1"].notna()
    ].copy()
    parsed = parse_flora_utm_1x1(flora_raw["UTM 1x1"])
    flora_xy = pd.concat([flora_raw.reset_index(drop=True), parsed], axis=1)
    flora_xy = flora_xy[flora_xy["UTMCODE1X1"].notna()].copy()

    lookup = (
        utm_net.drop(columns="geometry", errors="ignore")[
            ["UTMCODE1X1", "XCENTROIDE", "YCENTROIDE"]
        ]
        .drop_duplicates("UTMCODE1X1")
    )
    flora_xy = flora_xy.merge(lookup, on="UTMCODE1X1", how="left")
    flora_xy["XCENTROIDE"] = flora_xy["XCENTROIDE"].where(
        flora_xy["XCENTROIDE"].notna(),
        pd.to_numeric(flora_xy["utm_e"], errors="coerce") * 1000 + 500,
    )
    flora_xy["YCENTROIDE"] = flora_xy["YCENTROIDE"].where(
        flora_xy["YCENTROIDE"].notna(),
        pd.to_numeric(flora_xy["utm_n"], errors="coerce") * 1000 + 500,
    )
    flora_xy = flora_xy[flora_xy["XCENTROIDE"].notna() & flora_xy["YCENTROIDE"].notna()]
    dx, dy = offset_500_xy(flora_xy["quad_500"])
    flora_xy["XCENTROIDE"] = flora_xy["XCENTROIDE"] + dx.values
    flora_xy["YCENTROIDE"] = flora_xy["YCENTROIDE"] + dy.values
    flora_xy["grid_m"] = flora_xy["quad_500"].notna().map({True: 500, False: 1000})
    flora_xy["coordinateUncertaintyInMeters"] = flora_xy["grid_m"].map(
        {500: UNC_500M, 1000: UNC_1KM}
    )
    flora_xy["class"] = flora_xy["UTMCODE1X1"].map(
        lambda c: "PAC" if c in pac_codes else "PAU"
    )
    flora_xy["source"] = "Biel"
    flora_xy["scientificName"] = flora_xy["Especie"]
    flora_xy["eventDate"] = pd.to_datetime(flora_xy["Data"], errors="coerce", utc=True)

    gdf = gpd.GeoDataFrame(
        flora_xy[
            [
                "eventDate",
                "scientificName",
                "class",
                "coordinateUncertaintyInMeters",
                "source",
            ]
        ],
        geometry=[
            Point(xy) for xy in zip(flora_xy["XCENTROIDE"], flora_xy["YCENTROIDE"])
        ],
        crs=CRS,
    )
    pac = gdf[gdf["class"] == "PAC"][MERGE_COLS].copy()
    pau = gdf[gdf["class"] == "PAU"][MERGE_COLS].copy()
    print(f"Biel Mallorca: PAC={len(pac)} PAU={len(pau)} pac_cells={len(pac_codes)}")
    return pac, pau


def gbif_get(url: str, retries: int = 5) -> dict | None:
    for i in range(retries):
        try:
            r = SESSION.get(url, timeout=60)
            if r.status_code == 429:
                time.sleep(2 ** i)
                continue
            r.raise_for_status()
            return r.json()
        except Exception as exc:  # noqa: BLE001
            if i == retries - 1:
                print(f"GBIF fail {url}: {exc}")
                return None
            time.sleep(1.5 * (i + 1))
    return None


def resolve_taxon(name: str) -> dict:
    empty = {
        "initial_name": name,
        "accepted_key": pd.NA,
        "accepted_rank": "OTHER",
        "accepted_name": pd.NA,
        "accepted_species_key": pd.NA,
        "accepted_species_name": pd.NA,
        "accepted_order": pd.NA,
        "accepted_family": pd.NA,
    }
    match = gbif_get(
        "https://api.gbif.org/v1/species/match?"
        f"name={quote(name)}&rank=SPECIES"
    )
    if not match or not match.get("usageKey"):
        return empty
    rank = match.get("rank") or "OTHER"
    if rank in {"GENUS", "FAMILY", "CLASS"}:
        empty["accepted_rank"] = rank
        empty["accepted_order"] = match.get("order")
        empty["accepted_family"] = match.get("family")
        return empty

    accepted_key = match.get("acceptedUsageKey") or match.get("usageKey")
    acc = gbif_get(f"https://api.gbif.org/v1/species/{accepted_key}")
    if not acc:
        return empty
    accepted_name = acc.get("scientificName")
    acc_rank = acc.get("rank")
    if acc_rank == "SPECIES":
        species_key = accepted_key
        species_name = accepted_name
    elif acc.get("speciesKey"):
        sp = gbif_get(f"https://api.gbif.org/v1/species/{acc['speciesKey']}")
        if not sp:
            species_key, species_name = pd.NA, pd.NA
        else:
            species_key = sp.get("key")
            species_name = sp.get("scientificName")
    else:
        species_key, species_name = pd.NA, pd.NA

    return {
        "initial_name": name,
        "accepted_key": accepted_key,
        "accepted_rank": rank,
        "accepted_name": accepted_name,
        "accepted_species_key": species_key,
        "accepted_species_name": species_name,
        "accepted_order": match.get("order"),
        "accepted_family": match.get("family"),
    }


def build_taxonkey(taxa: list[str], review_path: Path) -> pd.DataFrame:
    old_review = pd.read_excel(review_path) if review_path.exists() else pd.DataFrame()
    by_name = {}
    if len(old_review):
        for _, row in old_review.iterrows():
            by_name[normalize_scientific_name(row["initial_name"])] = row.to_dict()

    rows = []
    missing = [t for t in taxa if t not in by_name]
    print(f"taxonkey: reuse={len(taxa) - len(missing)} fetch={len(missing)} total={len(taxa)}")
    for i, name in enumerate(taxa, 1):
        if name in by_name:
            rows.append(by_name[name])
            continue
        if i % 50 == 0 or i == 1:
            print(f"  GBIF {i}/{len(taxa)} ({name})")
        rows.append(resolve_taxon(name))
        time.sleep(0.05)
    out = pd.DataFrame(rows)
    # stable column order
    cols = [
        "initial_name",
        "accepted_key",
        "accepted_rank",
        "accepted_name",
        "accepted_species_key",
        "accepted_species_name",
        "accepted_order",
        "accepted_family",
    ]
    for c in cols:
        if c not in out.columns:
            out[c] = pd.NA
    return out[cols]


def clean_taxonkey(taxonkey: pd.DataFrame) -> pd.DataFrame:
    bad_orders = {
        "Bryales", "Bryopsidales", "Ceramiales", "Cladophorales",
        "Corallinales", "Dasycladales", "Dicranales", "Fossombroniales", "Funariales",
        "Hypnales", "Jungermanniales", "Marchantiales", "Nemaliales", "Pelliales",
        "Peyssonneliales", "Porellales", "Pottiales", "Scouleriales", "Siphonocladales",
    }
    bad_families = {
        "Posidoniaceae", "Potamogetonaceae", "Ulvaceae", "Phyllophoraceae",
        "Ceratophyllaceae", "Anadyomenaceae", "Zosteraceae", "Ruppiaceae",
        "Bonnemaisoniaceae", "Cymodoceaceae", "Hydrocharitaceae",
    }
    tk = taxonkey.copy()
    tk = tk[~tk["accepted_rank"].isin(["CLASS", "OTHER", "GENUS", "FAMILY"])]
    tk = tk[tk["accepted_species_key"].notna()]
    tk = tk[~tk["accepted_order"].isin(bad_orders)]
    tk = tk[~tk["accepted_family"].isin(bad_families)]
    tk = tk[tk["accepted_key"].notna()]
    return tk


def main() -> None:
    print("Loading Flora_net…")
    utm_net = gpd.read_file(RAW / "Flora_net.gpkg")

    print("Loading sources…")
    biodibal = load_biodibal()
    flora_pac, flora_pau, gbif_po = load_gbif_and_flora(utm_net)
    biel_pac, biel_pau = load_biel(utm_net)

    flora_pac = pd.concat([flora_pac, biel_pac], ignore_index=True)
    flora_pac = gpd.GeoDataFrame(flora_pac, geometry="geometry", crs=CRS)
    flora_pau = pd.concat([flora_pau, biel_pau], ignore_index=True)
    flora_pau = gpd.GeoDataFrame(flora_pau, geometry="geometry", crs=CRS)

    full_data = pd.concat([flora_pac, flora_pau, biodibal, gbif_po], ignore_index=True)
    full_data = gpd.GeoDataFrame(full_data, geometry="geometry", crs=CRS)
    full_data["scientificName_original"] = full_data["scientificName"]
    full_data["scientificName"] = full_data["scientificName"].map(normalize_scientific_name)
    print("counts by class/source:")
    print(full_data.groupby(["class", "source"], dropna=False).size())

    taxa = sorted(full_data["scientificName"].dropna().unique())
    review_path = RAW / "taxonkey_list_review.xlsx"
    taxonkey = build_taxonkey(taxa, review_path)
    taxonkey_path = RAW / "taxonkey_list.xlsx"
    taxonkey.to_excel(taxonkey_path, index=False)
    print(f"wrote {taxonkey_path} rows={len(taxonkey)}")

    # extend review with any new initial_name (keep prior manual rows)
    if review_path.exists():
        old = pd.read_excel(review_path)
        old_names = set(old["initial_name"].map(normalize_scientific_name))
        new_rows = taxonkey[~taxonkey["initial_name"].isin(old_names)]
        review_out = pd.concat([old, new_rows], ignore_index=True)
    else:
        review_out = taxonkey
    review_out.to_excel(review_path, index=False)
    print(f"wrote {review_path} rows={len(review_out)}")

    tk_clean = clean_taxonkey(taxonkey)
    name_lookup = tk_clean[
        [
            "initial_name",
            "accepted_key",
            "accepted_name",
            "accepted_species_key",
            "accepted_species_name",
        ]
    ]
    taxa_accepted = (
        tk_clean.groupby("accepted_species_key", as_index=False)
        .agg(accepted_name_canonical=("accepted_species_name", "first"))
    )

    merged = full_data.merge(
        name_lookup, left_on="scientificName", right_on="initial_name", how="left"
    )
    merged = merged[merged["accepted_species_key"].notna()].copy()
    merged = merged.merge(taxa_accepted, on="accepted_species_key", how="left")
    merged["scientificName"] = merged["accepted_name_canonical"]
    out = merged[
        [
            "eventDate",
            "scientificName",
            "class",
            "coordinateUncertaintyInMeters",
            "source",
            "scientificName_original",
            "geometry",
        ]
    ].copy()
    out = gpd.GeoDataFrame(out, geometry="geometry", crs=CRS)

    out_path = RAW / "full_data_saez.gpkg"
    if out_path.exists():
        out_path.unlink()
    out.to_file(out_path, driver="GPKG")
    print(
        f"wrote {out_path} rows={len(out)} species={out['scientificName'].nunique()} "
        f"sources={out['source'].value_counts().to_dict()}"
    )


if __name__ == "__main__":
    main()
