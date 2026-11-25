#!/usr/bin/env python3
"""Cross dataset taxa vs Clave (Mallorca) and Alien Flora Balearic list.

Output columns (main):
  en_clave / nombre_en_clave
  en_aloctonas / nombre_en_aloctonas
  eliminar_candidato = SI when neither source has the taxon

Also sheet rescate_texto: for eliminar_candidatos, search Quarta text,
clave_all (other islands), original-name variants, and fuzzy alien.
Does NOT auto-flip en_clave/en_aloctonas — only proposes for review.

Clave joins reuse matches from clave_mallorca_vs_dataset.xlsx.
Alien joins by canonical binomial (exact, then species-level collapse).
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
from rapidfuzz import fuzz, process

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from extract_clave_vs_dataset import (  # noqa: E402
    canonical_binomial,
    parse_genus_header,
    species_level_key,
    title_case_genus,
)

DATASET = ROOT / "data" / "raw" / "species_list_review.xlsx"
CLAVE_MATCH = ROOT / "data" / "raw" / "clave_mallorca_vs_dataset.xlsx"
CLAVE_TAXA = ROOT / "data" / "raw" / "clave_mallorca_especies.xlsx"
QUARTA_TXT = ROOT / "data" / "raw" / "Quarta_edicio.txt"
ALIEN = ROOT / "data" / "raw" / "Alien_Flora_Balearic_Records.xlsx"
OUT = ROOT / "data" / "raw" / "dataset_clave_aloctonas_review.xlsx"

ISLAND_NEAR_RE = re.compile(
    r"\(([^)]*\b(?:Mll|Mn|Ei|Fo)\b[^)]*)\)", re.IGNORECASE
)

SOLID_CLAVE = {
    "gbif_species_key",
    "name_exact",
    "name_exact_species_level",
    "name_exact_no_gbif",
    "name_exact_gbif_miss",
    "clave_paren_synonym",
    "fuzzy_name_fallback",
}


def load_clave_hooks(matches: pd.DataFrame) -> dict[int, dict]:
    m_ok = matches[
        matches["match_type"].isin(SOLID_CLAVE)
        & matches["dataset_accepted_species_key"].notna()
    ].copy()
    m_ok["prio"] = m_ok["match_type"].map(
        lambda t: 0 if t == "fuzzy_name_fallback" else 1
    )
    m_ok = m_ok.sort_values(["prio", "match_score"], ascending=[False, False])
    by_key: dict[int, dict] = {}
    for _, r in m_ok.iterrows():
        k = int(r["dataset_accepted_species_key"])
        if k in by_key:
            continue
        by_key[k] = {
            "en_clave": "SI",
            "nombre_en_clave": r["taxon_clave"],
            "clave_match_type": r["match_type"],
            "clave_match_note": r.get("match_note", "") or "",
            "clave_dubious": bool(r.get("dubious_hook", False)),
        }
    return by_key


def load_alien(path: Path) -> tuple[dict, dict]:
    alien = pd.read_excel(path, sheet_name="Sheet1")
    alien.columns = [c.strip() if isinstance(c, str) else c for c in alien.columns]
    alien = alien.copy()
    alien["Taxa"] = (
        alien["Taxa"]
        .astype(str)
        .str.replace("\xa0", " ", regex=False)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )
    alien["canonical"] = alien["Taxa"].map(canonical_binomial)
    alien["species_key_name"] = alien["canonical"].map(species_level_key)
    alien["Ma"] = pd.to_numeric(alien.get("Ma"), errors="coerce")

    by_canon: dict = {}
    by_sp: dict = defaultdict(list)
    for _, r in alien.iterrows():
        c = r["canonical"]
        if not c or (isinstance(c, float) and pd.isna(c)):
            continue
        if c not in by_canon:
            by_canon[c] = r
        sk = r["species_key_name"]
        if sk:
            by_sp[sk].append(r)
    return by_canon, by_sp


def match_alien(canon, spkey, by_canon, by_sp) -> dict:
    empty = {
        "en_aloctonas": "NO",
        "nombre_en_aloctonas": "",
        "aloctonas_match_type": "",
        "aloctonas_Ma": "",
        "aloctonas_family": "",
        "aloctonas_invasiveness": "",
    }
    hit = None
    mtype = ""
    if canon and canon in by_canon:
        hit = by_canon[canon]
        mtype = "name_exact"
    elif spkey and spkey in by_sp:
        rows = by_sp[spkey]
        exact = [x for x in rows if x["canonical"] == spkey]
        hit = (exact or rows)[0]
        mtype = "species_level"
    if hit is None:
        return empty
    ma = hit.get("Ma")
    return {
        "en_aloctonas": "SI",
        "nombre_en_aloctonas": hit["Taxa"],
        "aloctonas_match_type": mtype,
        "aloctonas_Ma": "SI" if pd.notna(ma) and float(ma) == 1 else "NO",
        "aloctonas_family": hit.get("Family", ""),
        "aloctonas_invasiveness": hit.get(
            "Invasiveness_Blackburn_simplified_analysis", ""
        ),
    }


def _genus_epi(canon: str | None) -> tuple[str, str] | None:
    if not canon or not isinstance(canon, str):
        return None
    parts = canon.split()
    if len(parts) < 2:
        return None
    return parts[0], parts[1]


def build_genus_line_index(text: str) -> list[tuple[int, str | None]]:
    """Map line start offsets -> current genus header context."""
    lines = text.splitlines(keepends=True)
    out: list[tuple[int, str | None]] = []
    current: str | None = None
    pos = 0
    for line in lines:
        raw = line.strip()
        head = parse_genus_header(raw) if raw else None
        if head and len(head) >= 2 and not head.endswith("ACEAE"):
            # skip obvious family-like / skip headers handled loosely
            if head.isupper() and len(head) >= 3:
                current = title_case_genus(head)
        out.append((pos, current))
        pos += len(line)
    return out


def genus_at(index: list[tuple[int, str | None]], char_pos: int) -> str | None:
    # binary-ish linear scan (text ~600k, lines ~few k — fine)
    cur = None
    for start, g in index:
        if start > char_pos:
            break
        cur = g
    return cur


def search_texto_clave(
    canon: str,
    text: str,
    genus_index: list[tuple[int, str | None]],
) -> dict:
    """Find Genus epithet or validated G. epithet under matching genus header."""
    ge = _genus_epi(canon)
    empty = {
        "texto_clave_hit": "NO",
        "texto_clave_form": "",
        "texto_clave_snippet": "",
        "texto_clave_islands_near": "",
        "texto_clave_genus_ctx": "",
    }
    if not ge:
        return empty
    genus, epi = ge
    patterns = [
        (rf"\b{re.escape(genus)}\s+(?:×\s*|x\s+)?{re.escape(epi)}\b", "full"),
        (rf"\b{re.escape(genus[0])}\.\s*(?:×\s*|x\s+)?{re.escape(epi)}\b", "abbrev1"),
    ]
    if len(genus) >= 2:
        patterns.append(
            (rf"\b{re.escape(genus[:2])}\.\s*(?:×\s*|x\s+)?{re.escape(epi)}\b", "abbrev2")
        )
    if len(genus) >= 3:
        patterns.append(
            (rf"\b{re.escape(genus[:3])}\.\s*(?:×\s*|x\s+)?{re.escape(epi)}\b", "abbrev3")
        )

    for pat, form in patterns:
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            ctx = genus_at(genus_index, m.start())
            # full genus name always ok; abbrev only if header genus matches
            if form != "full":
                if not ctx or ctx.lower() != genus.lower():
                    continue
            start = max(0, m.start() - 50)
            end = min(len(text), m.end() + 60)
            snip = re.sub(r"\s+", " ", text[start:end]).strip()
            islands = ""
            near = ISLAND_NEAR_RE.search(text[m.start() : m.end() + 80])
            if near:
                islands = near.group(1)
            return {
                "texto_clave_hit": "SI",
                "texto_clave_form": form,
                "texto_clave_snippet": snip[:180],
                "texto_clave_islands_near": islands,
                "texto_clave_genus_ctx": ctx or "",
            }
    return empty


def build_rescate_texto(elim: pd.DataFrame) -> pd.DataFrame:
    """Rescue search for taxa in neither clave-Mll join nor alien exact join."""
    if elim.empty:
        return elim.copy()

    clave_all = pd.read_excel(CLAVE_TAXA, sheet_name="clave_all")
    by_canon = {
        c: r
        for c, r in zip(clave_all["canonical"], clave_all.to_dict("records"))
        if isinstance(c, str)
    }
    by_sp: dict[str, list] = defaultdict(list)
    for c, r in by_canon.items():
        sk = species_level_key(c)
        if sk:
            by_sp[sk].append(r)
    clave_canons = list(by_canon.keys())
    by_genus_clave: dict[str, list[str]] = defaultdict(list)
    for c in clave_canons:
        by_genus_clave[c.split()[0].lower()].append(c)

    alien = pd.read_excel(ALIEN, sheet_name="Sheet1")
    alien.columns = [c.strip() if isinstance(c, str) else c for c in alien.columns]
    alien_canons: dict[str, str] = {}
    for t in alien["Taxa"].astype(str):
        t2 = re.sub(r"\s+", " ", t.replace("\xa0", " ")).strip()
        c = canonical_binomial(t2)
        if c and c not in alien_canons:
            alien_canons[c] = t2
    by_genus_alien: dict[str, list[str]] = defaultdict(list)
    for c in alien_canons:
        by_genus_alien[c.split()[0].lower()].append(c)

    variants = pd.read_excel(DATASET, sheet_name="gpkg_original_variants")
    orig_by_acc: dict[str, list[str]] = defaultdict(list)
    for _, r in variants.iterrows():
        acc = r["accepted_species_name"]
        orig = r["original_name"]
        if pd.isna(acc) or pd.isna(orig):
            continue
        if str(orig) != str(acc):
            orig_by_acc[str(acc)].append(str(orig))

    text = QUARTA_TXT.read_text(encoding="utf-8", errors="replace")
    genus_index = build_genus_line_index(text)

    rows = []
    for _, r in elim.iterrows():
        canon = r.get("canonical")
        sp = species_level_key(canon) if isinstance(canon, str) else None
        name = r["accepted_species_name"]

        # 1) exact / species-level in clave_all (often other islands)
        clave_hit = None
        clave_how = ""
        if isinstance(canon, str) and canon in by_canon:
            clave_hit = by_canon[canon]
            clave_how = "canonical_exact"
        elif sp and sp in by_sp:
            exact = [x for x in by_sp[sp] if x.get("canonical") == sp]
            clave_hit = (exact or by_sp[sp])[0]
            clave_how = "species_level"

        # 2) fuzzy same-genus vs clave_all
        fuzzy_clave_name = ""
        fuzzy_clave_score = None
        if isinstance(canon, str):
            pool = by_genus_clave.get(canon.split()[0].lower(), [])
            if pool:
                best = process.extractOne(canon, pool, scorer=fuzz.token_sort_ratio)
                if best and best[1] >= 90 and best[0] != canon:
                    # skip if already exact
                    if not clave_hit or best[0] != clave_hit.get("canonical"):
                        fuzzy_clave_name = best[0]
                        fuzzy_clave_score = float(best[1])

        # 3) text search validated
        texto = (
            search_texto_clave(canon, text, genus_index)
            if isinstance(canon, str)
            else {
                "texto_clave_hit": "NO",
                "texto_clave_form": "",
                "texto_clave_snippet": "",
                "texto_clave_islands_near": "",
                "texto_clave_genus_ctx": "",
            }
        )

        # 4) original-name variants → clave_all / text / alien
        via_orig = ""
        via_orig_target = ""
        via_orig_where = ""
        for orig in orig_by_acc.get(str(name), []):
            oc = canonical_binomial(orig)
            if not oc:
                continue
            if oc in by_canon:
                via_orig = orig
                via_orig_target = by_canon[oc].get("taxon_clave") or oc
                via_orig_where = "clave_all_via_original"
                break
            otext = search_texto_clave(oc, text, genus_index)
            if otext["texto_clave_hit"] == "SI":
                via_orig = orig
                via_orig_target = oc
                via_orig_where = "texto_clave_via_original"
                break
            if oc in alien_canons:
                via_orig = orig
                via_orig_target = alien_canons[oc]
                via_orig_where = "aloctonas_via_original"
                break

        # 5) fuzzy same-genus alien (high bar; typo rescue)
        fuzzy_alien_name = ""
        fuzzy_alien_score = None
        if isinstance(canon, str) and canon not in alien_canons:
            pool = by_genus_alien.get(canon.split()[0].lower(), [])
            if pool:
                best = process.extractOne(canon, pool, scorer=fuzz.token_sort_ratio)
                # require very high: crus-gallii vs crus-galli
                if best and best[1] >= 95:
                    fuzzy_alien_name = alien_canons[best[0]]
                    fuzzy_alien_score = float(best[1])

        en_clave_all = "SI" if clave_hit else "NO"
        islands = ""
        in_mll = ""
        nombre_clave_all = ""
        if clave_hit:
            islands = str(clave_hit.get("islands") or "")
            in_mll = "SI" if clave_hit.get("in_mallorca") else "NO"
            nombre_clave_all = str(
                clave_hit.get("taxon_clave") or clave_hit.get("canonical") or ""
            )

        # proposal for human
        if en_clave_all == "SI" and in_mll == "NO":
            propuesta = "en_clave_otras_islas"
        elif en_clave_all == "SI" and in_mll == "SI":
            propuesta = "rescate_clave_Mll_revisar"
        elif texto["texto_clave_hit"] == "SI":
            propuesta = "rescate_texto_clave_revisar"
        elif via_orig_where.startswith("clave") or via_orig_where.startswith("texto"):
            propuesta = "rescate_via_nombre_original_revisar"
        elif via_orig_where.startswith("aloctonas") or fuzzy_alien_name:
            propuesta = "rescate_aloctonas_revisar"
        elif fuzzy_clave_name:
            propuesta = "fuzzy_clave_revisar"
        else:
            propuesta = "sin_rescate_probable_error"

        rows.append(
            {
                "accepted_species_name": name,
                "canonical": canon,
                "accepted_species_key": r.get("accepted_species_key"),
                "n_obs_gpkg": r.get("n_obs_gpkg"),
                "PAC": r.get("PAC"),
                "PAU": r.get("PAU"),
                "POV": r.get("POV"),
                "POU": r.get("POU"),
                "propuesta": propuesta,
                "clave_all_hit": en_clave_all,
                "clave_all_how": clave_how,
                "nombre_en_clave_all": nombre_clave_all,
                "clave_all_islands": islands,
                "clave_all_in_mallorca": in_mll,
                "fuzzy_clave_nombre": fuzzy_clave_name,
                "fuzzy_clave_score": fuzzy_clave_score,
                **texto,
                "via_original_name": via_orig,
                "via_original_target": via_orig_target,
                "via_original_where": via_orig_where,
                "fuzzy_aloctonas_nombre": fuzzy_alien_name,
                "fuzzy_aloctonas_score": fuzzy_alien_score,
                "correccion_manual": "",
                "comentario": "",
            }
        )

    out = pd.DataFrame(rows)
    order = {
        "rescate_clave_Mll_revisar": 0,
        "rescate_texto_clave_revisar": 1,
        "rescate_via_nombre_original_revisar": 2,
        "en_clave_otras_islas": 3,
        "rescate_aloctonas_revisar": 4,
        "fuzzy_clave_revisar": 5,
        "sin_rescate_probable_error": 6,
    }
    out["_o"] = out["propuesta"].map(lambda x: order.get(x, 9))
    out = out.sort_values(["_o", "accepted_species_name"]).drop(columns="_o")
    return out


def main() -> None:
    if not ALIEN.exists():
        raise FileNotFoundError(f"Missing alien list: {ALIEN}")
    if not CLAVE_MATCH.exists():
        raise FileNotFoundError(
            f"Missing clave matches: {CLAVE_MATCH}. Run extract_clave_vs_dataset.py first."
        )
    if not CLAVE_TAXA.exists():
        raise FileNotFoundError(f"Missing clave taxa: {CLAVE_TAXA}")
    if not QUARTA_TXT.exists():
        raise FileNotFoundError(f"Missing Quarta text: {QUARTA_TXT}")

    acc = pd.read_excel(DATASET, sheet_name="gpkg_accepted").copy()
    acc["canonical"] = acc["accepted_species_name"].map(canonical_binomial)
    acc["species_key_name"] = acc["canonical"].map(species_level_key)
    acc["accepted_species_key"] = pd.to_numeric(
        acc["accepted_species_key"], errors="coerce"
    ).astype("Int64")

    matches = pd.read_excel(CLAVE_MATCH, sheet_name="matches")
    clave_by_key = load_clave_hooks(matches)
    alien_by_canon, alien_by_sp = load_alien(ALIEN)

    rows = []
    for _, r in acc.iterrows():
        k = r["accepted_species_key"]
        if pd.notna(k) and int(k) in clave_by_key:
            cinfo = clave_by_key[int(k)]
        else:
            cinfo = {
                "en_clave": "NO",
                "nombre_en_clave": "",
                "clave_match_type": "",
                "clave_match_note": "",
                "clave_dubious": False,
            }
        ainfo = match_alien(
            r["canonical"], r["species_key_name"], alien_by_canon, alien_by_sp
        )
        en_alguno = (
            "SI"
            if cinfo["en_clave"] == "SI" or ainfo["en_aloctonas"] == "SI"
            else "NO"
        )
        rows.append(
            {
                "accepted_species_name": r["accepted_species_name"],
                "canonical": r["canonical"],
                "accepted_species_key": k,
                "n_obs_gpkg": r["n_obs_gpkg"],
                "PAC": r["PAC"],
                "PAU": r["PAU"],
                "POV": r["POV"],
                "POU": r["POU"],
                "en_clave": cinfo["en_clave"],
                "nombre_en_clave": cinfo["nombre_en_clave"],
                "en_aloctonas": ainfo["en_aloctonas"],
                "nombre_en_aloctonas": ainfo["nombre_en_aloctonas"],
                "en_clave_o_aloctonas": en_alguno,
                "eliminar_candidato": "SI" if en_alguno == "NO" else "NO",
                "clave_match_type": cinfo["clave_match_type"],
                "clave_dubious": "SI" if cinfo["clave_dubious"] else "NO",
                "clave_match_note": cinfo["clave_match_note"],
                "aloctonas_match_type": ainfo["aloctonas_match_type"],
                "aloctonas_registrada_Mallorca": ainfo["aloctonas_Ma"],
                "aloctonas_family": ainfo["aloctonas_family"],
                "aloctonas_invasiveness": ainfo["aloctonas_invasiveness"],
                "correccion_manual_en_clave": "",
                "correccion_manual_en_aloctonas": "",
                "comentario": "",
            }
        )

    out = pd.DataFrame(rows)
    out["_ord"] = out["eliminar_candidato"].map({"SI": 0, "NO": 1})
    out = out.sort_values(
        ["_ord", "en_clave", "en_aloctonas", "accepted_species_name"]
    ).drop(columns="_ord")

    summary = pd.DataFrame(
        [
            {"metric": "dataset_taxa", "value": len(out)},
            {"metric": "en_clave_SI", "value": int((out["en_clave"] == "SI").sum())},
            {
                "metric": "en_aloctonas_SI",
                "value": int((out["en_aloctonas"] == "SI").sum()),
            },
            {
                "metric": "en_ambos",
                "value": int(
                    ((out["en_clave"] == "SI") & (out["en_aloctonas"] == "SI")).sum()
                ),
            },
            {
                "metric": "solo_clave",
                "value": int(
                    ((out["en_clave"] == "SI") & (out["en_aloctonas"] == "NO")).sum()
                ),
            },
            {
                "metric": "solo_aloctonas",
                "value": int(
                    ((out["en_clave"] == "NO") & (out["en_aloctonas"] == "SI")).sum()
                ),
            },
            {
                "metric": "en_ninguno_eliminar",
                "value": int((out["eliminar_candidato"] == "SI").sum()),
            },
            {
                "metric": "clave_dubious",
                "value": int((out["clave_dubious"] == "SI").sum()),
            },
            {"metric": "source_clave_matches", "value": CLAVE_MATCH.name},
            {"metric": "source_alien", "value": ALIEN.name},
            {
                "metric": "source_dataset",
                "value": "species_list_review.xlsx#gpkg_accepted",
            },
        ]
    )

    elim = out[out["eliminar_candidato"] == "SI"].copy()
    solo_alien = out[
        (out["en_clave"] == "NO") & (out["en_aloctonas"] == "SI")
    ].copy()
    print(f"Building rescate_texto for {len(elim)} eliminar_candidatos...")
    rescate = build_rescate_texto(elim)
    rescate_counts = (
        rescate["propuesta"].value_counts().rename_axis("metric").reset_index(name="value")
        if not rescate.empty
        else pd.DataFrame(columns=["metric", "value"])
    )
    rescate_counts["metric"] = "rescate_" + rescate_counts["metric"].astype(str)
    summary = pd.concat([summary, rescate_counts], ignore_index=True)

    compact = out[
        [
            "accepted_species_name",
            "en_clave",
            "nombre_en_clave",
            "en_aloctonas",
            "nombre_en_aloctonas",
            "eliminar_candidato",
            "n_obs_gpkg",
            "PAC",
            "PAU",
            "POV",
            "POU",
            "correccion_manual_en_clave",
            "correccion_manual_en_aloctonas",
            "comentario",
        ]
    ]

    with pd.ExcelWriter(OUT, engine="openpyxl") as xl:
        summary.to_excel(xl, sheet_name="summary", index=False)
        compact.to_excel(xl, sheet_name="para_corregir", index=False)
        out.to_excel(xl, sheet_name="revision", index=False)
        elim.to_excel(xl, sheet_name="eliminar_candidatos", index=False)
        solo_alien.to_excel(xl, sheet_name="solo_aloctonas", index=False)
        rescate.to_excel(xl, sheet_name="rescate_texto", index=False)

    print(summary.to_string(index=False))
    if not rescate.empty:
        print("\nrescate propuesta:")
        print(rescate["propuesta"].value_counts().to_string())
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
