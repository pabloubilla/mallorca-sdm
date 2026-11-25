#!/usr/bin/env python3
"""Extract Mallorca taxa from Clave.pdf (Gil & Llorens 2021) and match to dataset.

Primary join: GBIF accepted_species_key (synonym-safe: Otanthus <-> Achillea, etc.).
Name matching kept only as diagnostic / fallback.
"""

from __future__ import annotations

import os
import re
import subprocess
import time
import unicodedata
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests
from rapidfuzz import fuzz, process

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "data" / "raw" / "Clave.pdf"
DATASET_XLSX = ROOT / "data" / "raw" / "species_list_review.xlsx"
TAXONKEY_XLSX = ROOT / "data" / "raw" / "taxonkey_list.xlsx"
TAXONKEY_REVIEW_XLSX = ROOT / "data" / "raw" / "taxonkey_list_review.xlsx"
GBIF_CACHE_XLSX = ROOT / "data" / "raw" / "clave_gbif_cache.xlsx"
OUT_XLSX = ROOT / "data" / "raw" / "clave_mallorca_vs_dataset.xlsx"
CLAVE_ONLY_XLSX = ROOT / "data" / "raw" / "clave_mallorca_especies.xlsx"

SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "mallorca-sdm/extract_clave_vs_dataset"})

# Stop before alphabetical index / glossary (PDF page numbers)
INDEX_PDF_PAGE_START = 209

SKIP_HEADERS = {
    "SUMARI",
    "INTRODUCCIO",
    "INTRODUCCIÓ",
    "CLAUS",
    "GLOSSARI",
    "RECORDAR",
    "DICOTYLEDONES",
    "MONOCOTYLEDONES",
    "MAGNOLIOPSIDA",
    "TRACHEOPHYTINA",
    "SPERMATOPHYTINA",
    "LYCOPODIOPHYTINA",
    "PTERIDOPHYTINA",
    "GNETOPSIDA",
    "PINOPSIDA",
}

ISLAND_RE = re.compile(
    r"\((?P<islands>[^)]*\b(?:Mll|Mn|Ei|Fo)\b[^)]*)\)",
    re.IGNORECASE,
)

# Explicit non-naturalized marker in the key (often without island codes).
NO_NAT_RE = re.compile(r"\(\s*no\s+naturalitzad[ae]s?\s*\)", re.IGNORECASE)

# Strict scientific-name chunk (must START at abbrev or Latin genus, not Catalan).
TAXON_ABBREV_RE = re.compile(
    r"(?P<taxon>"
    # Require '.' so Catalan / "A contrallum" never match as abbrev.
    # "T.marklundii" (no space after dot) allowed via \s*.
    r"(?P<abbr>[A-ZÁÉÍÓÚÀÈÒÜÑÏ][a-záéíóúàèòüñïç]{0,4}\.)"
    r"\s*"
    r"(?:gr\.\s+|aggr\.?\s+)?"
    r"(?:×\s*)?[a-záéíóúàèòüñïç][a-záéíóúàèòüñïç\-]*"
    r"(?:\s+(?:subsp\.|ssp\.|var\.|f\.)\s+[a-záéíóúàèòüñïç][a-záéíóúàèòüñïç\-]*)?"
    r"(?:"
    r"\s*\([^)=*][^)]*\)"
    r"|\s*=\s*[^=(]+"
    r"|\s+(?:×\s*)?[A-Za-záéíóúàèòüñïç\.\-&']+\*?"
    r")*"
    r")\*?"
)

# Typo form without dot: "K cirrhosa" (only used as fallback for current genus).
TAXON_NODOT_ABBREV_RE = re.compile(
    r"(?P<taxon>"
    r"(?P<abbr>[A-ZÁÉÍÓÚÀÈÒÜÑÏ])\s+"
    r"(?:×\s*)?(?P<epi>[a-záéíóúàèòüñïç][a-záéíóúàèòüñïç\-]*)"
    r"(?:\s+(?:subsp\.|ssp\.|var\.|f\.)\s+[a-záéíóúàèòüñïç][a-záéíóúàèòüñïç\-]*)?"
    r"(?:"
    r"\s*\([^)=*][^)]*\)"
    r"|\s+(?:×\s*)?[A-Za-záéíóúàèòüñïç\.\-&']+\*?"
    r")*"
    r")\*?"
)

TAXON_FULLGENUS_RE = re.compile(
    r"(?P<taxon>"
    r"(?P<genus>[A-ZÁÉÍÓÚÀÈÒÜÑÏ][a-záéíóúàèòüñïç\-]+)"
    r"\s+(?:gr\.\s+|aggr\.?\s+)?"
    r"(?:×\s*)?[a-záéíóúàèòüñïç][a-záéíóúàèòüñïç\-]*"
    r"(?:\s+(?:subsp\.|ssp\.|var\.|f\.)\s+[a-záéíóúàèòüñïç][a-záéíóúàèòüñïç\-]*)?"
    r"(?:"
    r"\s*\([^)=*][^)]*\)"
    r"|\s*=\s*[^=(]+"
    r"|\s+(?:×\s*)?[A-Za-záéíóúàèòüñïç\.\-&']+\*?"
    r")*"
    r")\*?"
)

# Plain genus/family header, or with synonym note: GLEBIONIS (= CHRYSANTHEMUM)
GENUS_HEADER_RE = re.compile(
    r"^(?P<head>[A-ZÁÉÍÓÚÀÈÒÜÑÏÇ][A-ZÁÉÍÓÚÀÈÒÜÑÏÇ×\-]{1,40})"
    r"(?:\s*\(\s*=\s*[^)]*\))?\s*$"
)

FAMILY_RE = re.compile(r".*(ACEAE|PHYTINA|OPSIDA|ALES)$")

# Traditional family names without -aceae (still not genera)
FAMILY_ALIAS_HEADERS = {
    "COMPOSITAE",
    "CRUCIFERAE",
    "UMBELLIFERAE",
    "LABIATAE",
    "GRAMINEAE",
    "GUTTIFERAE",
    "LEGUMINOSAE",
}

CATALAN_FALSE_GENUS = {
    "ala", "ales", "altres", "amb", "angle", "angles", "antela", "anteres",
    "anvers", "apex", "apendixs", "aquenis", "base", "beines", "bractees",
    "branques", "capítol", "capitol", "cara", "cares", "color", "corol",
    "dors", "espata", "espores", "estams", "flors", "frondes", "fruit",
    "fulles", "gola", "inflorescencia", "inflorescencies", "label", "lamina",
    "marges", "ovaris", "pecíol", "peciol", "pecíols", "peciols", "pètals",
    "petals", "pinnes", "planta", "plantes", "presencia", "presència",
    "raquis", "segment", "segments", "sèpals", "sepals", "tiges", "tub",
}


def strip_accents(s: str) -> str:
    return "".join(
        c for c in unicodedata.normalize("NFKD", s) if not unicodedata.combining(c)
    )


def pdf_to_text(pdf: Path) -> str:
    return subprocess.check_output(
        ["pdftotext", "-layout", "-f", "1", "-l", str(INDEX_PDF_PAGE_START - 1), str(pdf), "-"],
        text=True,
        errors="replace",
    )


def title_case_genus(genus: str) -> str:
    genus = genus.strip()
    if not genus:
        return genus
    return genus[0].upper() + genus[1:].lower()


def parse_genus_header(raw: str) -> str | None:
    """Return header stem if line is a genus/family header (optional (= syn))."""
    m = GENUS_HEADER_RE.fullmatch(raw or "")
    return m.group("head") if m else None


def is_family_or_skip_header(head: str) -> bool:
    return bool(
        FAMILY_RE.match(head)
        or head in SKIP_HEADERS
        or strip_accents(head) in SKIP_HEADERS
        or head in FAMILY_ALIAS_HEADERS
    )


def genus_abbrevs(genus: str) -> set[str]:
    """Possible abbreviations for a genus header (A., Ch., Ph., ...)."""
    g = title_case_genus(genus)
    out = {f"{g[0]}."}
    # multi-letter prefixes used in the key (Ch., Ph., Oe., ...)
    for n in range(2, min(4, len(g)) + 1):
        out.add(f"{g[:n]}.")
    return out


def build_abbrev_index(genera: set[str]) -> dict[str, list[str]]:
    """Map 'R.' / 'Ch.' -> possible genera."""
    idx: dict[str, list[str]] = defaultdict(list)
    for g in genera:
        for ab in genus_abbrevs(g):
            idx[ab].append(g)
    # longer genera first when ambiguous
    for ab, gs in idx.items():
        idx[ab] = sorted(gs, key=len, reverse=True)
    return idx


def expand_taxon(
    raw: str,
    current_genus: str | None,
    genera: set[str],
    abbrev_index: dict[str, list[str]] | None = None,
) -> str:
    t = re.sub(r"\s+", " ", raw).strip(" .;,:")
    t = t.replace("× ", "×").replace(" ×", "×").rstrip("*")
    first = t.split()[0] if t else ""
    rest = t[len(first) :].lstrip() if first else ""

    # "K." or bare "K" abbreviation
    ab = first if first.endswith(".") else (first + "." if len(first) <= 5 and first.isalpha() else None)
    if ab and re.fullmatch(r"[A-Za-záéíóúàèòüñïç]{1,5}\.", ab):
        if current_genus and ab in genus_abbrevs(current_genus) and rest:
            return f"{title_case_genus(current_genus)} {rest}"
        if abbrev_index and ab in abbrev_index:
            gs = abbrev_index[ab]
            if current_genus and current_genus in gs:
                g = current_genus
            elif len(gs) == 1:
                g = gs[0]
            else:
                g = gs[0]
            if rest:
                return f"{title_case_genus(g)} {rest}"
        return t

    first_norm = first[0].upper() + first[1:].lower() if first else ""
    if first_norm in genera:
        return first_norm + t[len(first) :]
    if current_genus and t.lower().startswith(current_genus.lower() + " "):
        return title_case_genus(current_genus) + t[len(current_genus) :]
    return t


def name_for_canonical(expanded: str) -> str:
    """Drop synonymy / keep primary name for canonicalization."""
    t = re.sub(r"\s*\(=\s*[^)]*\)\s*", " ", expanded)
    t = re.sub(r"\s*\(\s*inclou\s+[^)]*\)\s*", " ", t, flags=re.IGNORECASE)
    t = re.sub(r"\s*\(\s*incl\.?\s+[^)]*\)\s*", " ", t, flags=re.IGNORECASE)
    t = re.sub(r"\s*=\s*.*$", "", t)
    # broken nested paren leftovers: "Font Quer) (inclou ..."
    t = re.sub(r"\)\s*$", "", t).strip()
    return re.sub(r"\s+", " ", t).strip(" .;,:")


def strip_synonym_notes(text: str) -> str:
    """Remove (= synonym) notes so pickers keep the primary name."""
    t = re.sub(r"\(\s*=\s*[^)]*\)", " ", text)
    return re.sub(r"\s+", " ", t).strip()


PAREN_SYN_RE = re.compile(r"\(\s*=\s*([^)]*)\)")


def parse_paren_synonyms(
    left: str,
    current_genus: str | None,
    genera: set[str],
    abbrev_index: dict[str, list[str]] | None = None,
) -> list[dict]:
    """Parse (= Synonym) notes into expanded/canonical synonym dicts."""
    out: list[dict] = []
    seen: set[str] = set()
    for m in PAREN_SYN_RE.finditer(left or ""):
        chunk = re.sub(r"\s+", " ", m.group(1)).strip(" .;,:")
        if not chunk:
            continue
        # Chunk is already the synonym body (no outer (=)).
        picked = pick_taxon_raw(chunk, current_genus, genera, abbrev_index)
        if not picked:
            # try as-is expand
            expanded = expand_taxon(chunk, current_genus, genera, abbrev_index)
        else:
            raw, hint = picked
            expanded = expand_taxon(raw, hint or current_genus, genera, abbrev_index)
        primary = name_for_canonical(expanded)
        canon = canonical_binomial(primary)
        if not canon:
            cleaned = re.sub(r"\([^)]*\)", " ", primary)
            canon = canonical_binomial(re.sub(r"\s+", " ", cleaned).strip())
        if not canon or canon in seen:
            continue
        seen.add(canon)
        out.append(
            {
                "synonym_raw": chunk,
                "synonym_expanded": expanded,
                "synonym_canonical": canon,
                "synonym_species_key": species_level_key(canon),
            }
        )
    return out

def extract_unlabeled_extras(text: str, genera: set[str], abbrev_index: dict[str, list[str]], by_canon: dict) -> list[dict]:
    """Capture hybrids / prose citations without standard (Mll,...) codes."""
    lines = text.splitlines()
    extras: list[dict] = []
    current_genus: str | None = None

    def add(expanded: str, islands: list[str], note: str, line_no: int, genus_ctx: str | None) -> None:
        primary = name_for_canonical(expanded)
        canon = canonical_binomial(primary)
        if not canon:
            cleaned = re.sub(r"\([^)]*\)", " ", primary)
            canon = canonical_binomial(re.sub(r"\s+", " ", cleaned).strip())
        if not canon or canon in by_canon:
            return
        g = canon.split()[0]
        # reject Catalan / non-genus garbage
        if g not in genera and g not in {"Helosciadium"} and g.lower() in CATALAN_FALSE_GENUS:
            return
        if g not in genera and g not in {"Helosciadium"}:
            # allow only if genus_ctx forced a known expansion
            if genus_ctx not in genera and genus_ctx not in {"Helosciadium"}:
                return
            if g[0].isupper() and g[1:].islower() and len(g) >= 4:
                pass  # Latin-looking forced genus (Helosciadium)
            else:
                return
        by_canon[canon] = {
            "taxon_clave_raw": expanded,
            "taxon_clave": expanded,
            "canonical": canon,
            "species_key": species_level_key(canon),
            "paren_synonyms": "",
            "paren_synonyms_raw": "",
            "rank_inferred": (
                "nothosubspecies"
                if "nothosubsp." in canon
                else "hybrid"
                if " x " in canon
                else "subspecies"
                if "subsp." in canon
                else "variety"
                if " var." in canon
                else "species"
            ),
            "islands": ",".join(x for x in ("Mll", "Mn", "Ei", "Fo") if x in islands),
            "in_mallorca": ("Mll" in islands) or (not islands and "Mallorca" in note),
            "no_naturalitzada": False,
            "genus_context": genus_ctx,
            "family_context": None,
            "source_line": line_no,
            "extraction_note": note,
        }
        extras.append(by_canon[canon])

    # Asplenium-style hybrid lists and inline hybrids
    hybrid_line = re.compile(
        r"^((?:[A-ZÁÉÍÓÚÀÈÒÜÑÏ][a-záéíóúàèòüñïç]{0,4}\.|[A-ZÁÉÍÓÚÀÈÒÜÑÏ][a-záéíóúàèòüñïç\-]+)"
        r"\s*(?:x|×)\s*[a-záéíóúàèòüñïç][a-záéíóúàèòüñïç\-]*"
        r"(?:\s+nothosubsp\.\s+[a-záéíóúàèòüñïç\-]+)?)"
        r"(?:\s|\(|$)",
        re.IGNORECASE,
    )
    # Species with nothosubsp. without leading x
    notho_line = re.compile(
        r"^((?:[A-ZÁÉÍÓÚÀÈÒÜÑÏ][a-záéíóúàèòüñïç]{0,4}\.|[A-ZÁÉÍÓÚÀÈÒÜÑÏ][a-záéíóúàèòüñïç\-]+)"
        r"\s+[a-záéíóúàèòüñïç][a-záéíóúàèòüñïç\-]*"
        r"\s+nothosubsp\.\s+[a-záéíóúàèòüñïç\-]+)",
        re.IGNORECASE,
    )

    for i, line in enumerate(lines):
        raw = line.strip()
        head = parse_genus_header(raw)
        if head:
            if not is_family_or_skip_header(head):
                current_genus = title_case_genus(head)
            continue
        if not raw:
            continue

        for rx in (hybrid_line, notho_line):
            m = rx.match(raw)
            if not m:
                continue
            if not current_genus:
                continue
            expanded = expand_taxon(m.group(1), current_genus, genera, abbrev_index)
            # require expansion into current genus
            if not expanded.lower().startswith(current_genus.lower()):
                continue
            note = "hybrid_list_no_island_code"
            islands = ["Mll"] if current_genus == "Asplenium" else []
            add(expanded, islands, note, i + 1, current_genus)

    joined = re.sub(r"\s+", " ", text)

    # Prose: "cita de S. multiceps Coss. & Durieu a Mallorca"
    for m in re.finditer(
        r"cita de\s+(S\.\s*multiceps\b(?:\s+[A-Za-záéíóúàèòüñïç\.\-&]+)*?)\s+a\s+Mallorca",
        joined,
        flags=re.IGNORECASE,
    ):
        expanded = expand_taxon(m.group(1).strip(), "Sedum", genera, abbrev_index)
        add(expanded, ["Mll"], "prose_citation_mallorca", 0, "Sedum")

    # Other named hybrids in prose with island context nearby
    # C. x sinaica (Crataegus), L. x strigocamara (Lantana)
    for m in re.finditer(
        r"\((C\.\s*(?:x|×)\s*sinaica\b[^)]{0,40})\)",
        joined,
        flags=re.IGNORECASE,
    ):
        expanded = expand_taxon(m.group(1).strip(), "Crataegus", genera, abbrev_index)
        add(expanded, ["Mll"], "prose_hybrid_crataegus", 0, "Crataegus")
    for m in re.finditer(
        r"(L\.\s*(?:x|×)\s*strigocamara\b(?:\s+[A-Za-záéíóúàèòüñïç\.\-&]+)*)",
        joined,
        flags=re.IGNORECASE,
    ):
        expanded = expand_taxon(m.group(1).strip(), "Lantana", genera, abbrev_index)
        add(expanded, ["Mll"], "prose_hybrid_lantana", 0, "Lantana")

    # Prose hybrid under Apium (= Helosciadium)
    for m in re.finditer(
        r"\((H\.\s*(?:x|×)\s*clandestinum\b[^)]{0,60})\)",
        joined,
        flags=re.IGNORECASE,
    ):
        # Name published as Helosciadium hybrid
        expanded = "Helosciadium x clandestinum"
        add(expanded, ["Mn"], "prose_hybrid_menorca", 0, "Helosciadium")

    return extras


def pick_taxon_raw(
    left: str,
    current_genus: str | None,
    genera: set[str],
    abbrev_index: dict[str, list[str]] | None = None,
) -> tuple[str, str | None] | None:
    """Pick rightmost plausible scientific name. Returns (raw, resolved_genus_hint)."""
    # Keep primary: "G. communis L. (= G. illyricus ...)" must not pick the synonym.
    left = strip_synonym_notes(left)
    ranked: list[tuple[int, int, str, str | None]] = []

    def _rank_abbrev(raw: str, ab: str, start: int) -> None:
        nonlocal ranked
        if not ab.endswith("."):
            ab = ab + "."
        if start > 0 and left[start - 1].isalpha():
            return
        hint = None
        score = 0
        if current_genus and ab in genus_abbrevs(current_genus):
            score = 3
            hint = current_genus
        elif abbrev_index and ab in abbrev_index:
            gs = abbrev_index[ab]
            if len(gs) == 1:
                score = 2
                hint = gs[0]
            elif current_genus and current_genus in gs:
                score = 2
                hint = current_genus
            else:
                score = 1
                hint = gs[0]
        else:
            return
        expanded = expand_taxon(raw, hint or current_genus, genera, abbrev_index)
        primary = name_for_canonical(expanded)
        if not canonical_binomial(primary):
            cleaned = re.sub(r"\([^)]*\)", " ", primary)
            cleaned = re.sub(r"\s+", " ", cleaned).strip()
            if not canonical_binomial(cleaned):
                return
        ranked.append((score, start, raw, hint))

    for m in TAXON_ABBREV_RE.finditer(left):
        _rank_abbrev(m.group("taxon").strip().rstrip("*"), m.group("abbr"), m.start("abbr"))

    # fallback: "K cirrhosa" only if matches current genus initial and no dotted hit yet
    if current_genus and not any(s >= 2 for s, _, _, _ in ranked):
        initial = current_genus[0]
        for m in TAXON_NODOT_ABBREV_RE.finditer(left):
            if m.group("abbr") != initial:
                continue
            # avoid Catalan "A contrallum" / "A la"
            epi = m.group("epi").lower()
            if epi in CATALAN_FALSE_GENUS or epi in {
                "contrallum", "la", "les", "el", "els", "un", "una", "amb", "del", "de",
            }:
                continue
            _rank_abbrev(m.group("taxon").strip().rstrip("*"), m.group("abbr"), m.start("abbr"))

    for m in TAXON_FULLGENUS_RE.finditer(left):
        raw = m.group("taxon").strip().rstrip("*")
        g = m.group("genus")
        g_norm = g[0].upper() + g[1:].lower()
        if g_norm.lower() in CATALAN_FALSE_GENUS:
            continue
        if g_norm == current_genus:
            score = 3
        elif g_norm in genera:
            score = 2
        elif (
            len(g_norm) >= 4
            and g_norm[0].isupper()
            and g_norm[1:].islower()
            and g_norm.lower() not in CATALAN_FALSE_GENUS
        ):
            # prose citation of a genus not used as section header
            score = 1
        else:
            continue
        expanded = expand_taxon(raw, g_norm, genera, abbrev_index)
        primary = name_for_canonical(expanded)
        if not canonical_binomial(primary):
            cleaned = re.sub(r"\([^)]*\)", " ", primary)
            cleaned = re.sub(r"\s+", " ", cleaned).strip()
            if not canonical_binomial(cleaned):
                continue
        ranked.append((score, m.start(), raw, g_norm))

    if not ranked:
        return None
    ranked.sort(key=lambda x: (x[0], x[1]))
    best = ranked[-1]
    return best[2], best[3]


def canonical_binomial(name: str) -> str | None:
    """Genus + epithet (+ optional infraspecific / hybrid mark) without authors."""
    if not name:
        return None
    s = strip_accents(name)
    s = s.replace("*", "")
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\baggr\.?\b", " ", s, flags=re.IGNORECASE)
    s = re.sub(r"\bgr\.\b", " ", s, flags=re.IGNORECASE)
    s = s.replace("×", "x ")
    s = re.sub(r"\s+", " ", s).strip()

    m = re.match(
        r"^([A-Za-z][A-Za-z\-]+)\s+(x\s+)?([a-z][a-z\-]+)",
        s,
    )
    if not m:
        return None
    genus, hybrid, epi = m.group(1), m.group(2), m.group(3)
    genus = genus[0].upper() + genus[1:].lower()
    base = f"{genus} {'x ' if hybrid else ''}{epi.lower()}"

    # nothosubsp. / subsp. / var. after authorship or directly
    infra = re.search(
        r"\b(nothosubsp\.|subsp\.|ssp\.|var\.|f\.)\s+([a-z][a-z\-]+)",
        s[m.end() :],
        flags=re.IGNORECASE,
    )
    if infra:
        rank = infra.group(1).lower()
        infra_epi = infra.group(2).lower()
        if rank.startswith("nothosub"):
            rank_n = "nothosubsp."
        elif rank.startswith("ss") or rank.startswith("sub"):
            rank_n = "subsp."
        elif rank.startswith("var"):
            rank_n = "var."
        else:
            rank_n = "f."
        base = f"{base} {rank_n} {infra_epi}"
    return base


def species_level_key(canon: str | None) -> str | None:
    """Collapse to genus + species epithet (drop infraspecific)."""
    if not canon:
        return None
    parts = canon.split()
    if len(parts) < 2:
        return None
    return f"{parts[0]} {parts[1]}"


def parse_islands(raw: str) -> list[str]:
    parts = re.split(r"[,;/]+", raw)
    out = []
    for p in parts:
        p = p.strip().rstrip(".")
        if not p:
            continue
        # strip uncertainty / extinct marks: Mll†, Mll?, Mn?
        pl = re.sub(r"[†?\s]+", "", p)
        low = pl.lower()
        if low.startswith("mll"):
            out.append("Mll")
        elif low.startswith("mn"):
            out.append("Mn")
        elif low.startswith("ei"):
            out.append("Ei")
        elif low.startswith("fo"):
            out.append("Fo")
        else:
            out.append(p)
    return out


def load_clave_text(source: Path | None = None) -> tuple[str, str]:
    """Load Quarta .doc/.txt preferentially; fall back to Clave.pdf.

    Returns (text, source_label).
    """
    quarta_txt = ROOT / "data" / "raw" / "Quarta_edicio.txt"
    quarta_doc = ROOT / "data" / "raw" / "Quarta edició.doc"
    antiword = ROOT / ".tools" / "antiword" / "antiword"
    mapping = ROOT / ".tools" / "antiword" / "UTF-8.txt"

    if source is not None:
        paths = [source]
    else:
        paths = [quarta_txt, quarta_doc, PDF]

    for path in paths:
        if not path.exists():
            continue
        if path.suffix.lower() == ".txt":
            text = path.read_text(encoding="utf-8", errors="replace")
            # drop glossary; keep body (TOC also mentions Índex — use last GLOSSARI)
            cut = text.rfind("GLOSSARI")
            if cut > 0:
                text = text[:cut]
            return text, path.name
        if path.suffix.lower() == ".doc" and antiword.exists():
            env = os.environ.copy()
            env["ANTIWORDHOME"] = str(antiword.parent)
            text = subprocess.check_output(
                [str(antiword), "-m", str(mapping), str(path)],
                text=True,
                errors="replace",
                env=env,
            )
            cut = text.rfind("GLOSSARI")
            if cut > 0:
                text = text[:cut]
            # cache txt for faster re-runs
            quarta_txt.write_text(text, encoding="utf-8")
            return text, path.name
        if path.suffix.lower() == ".pdf":
            return pdf_to_text(path), path.name
    raise FileNotFoundError("No Clave source found (Quarta doc/txt or Clave.pdf)")


def collect_genera(lines: list[str]) -> set[str]:
    genera: set[str] = set()
    for line in lines:
        head = parse_genus_header(line.strip())
        if not head or is_family_or_skip_header(head):
            continue
        genera.add(title_case_genus(head))
    return genera


def extract_clave_taxa(text: str, mallorca_only: bool = False) -> pd.DataFrame:
    lines = text.splitlines()
    genera = collect_genera(lines)
    abbrev_index = build_abbrev_index(genera)
    current_genus: str | None = None
    current_family: str | None = None
    buf: list[str] = []
    records: list[dict] = []
    # dedupe by canonical; merge island sets
    by_canon: dict[str, dict] = {}

    for i, line in enumerate(lines):
        raw = line.strip()
        # Keep buffer across blank/whitespace lines: antiword often inserts
        # spacer lines between a split authority and the island codes.
        if not raw:
            continue

        if re.fullmatch(r"\d{1,3}", raw):
            continue

        head = parse_genus_header(raw)
        if head:
            if is_family_or_skip_header(head):
                if head.endswith("ACEAE") or head in FAMILY_ALIAS_HEADERS:
                    current_family = title_case_genus(head)
                buf = []
                continue
            current_genus = title_case_genus(head)
            buf = []
            continue

        buf.append(raw)
        if len(buf) > 12:
            buf = buf[-12:]
        window = " ".join(buf[-6:])
        window = re.sub(r"\.{2,}", " ", window)
        window = re.sub(r"\s+", " ", window).strip()

        def _record_from_left(
            left: str,
            islands: list[str],
            *,
            no_naturalitzada: bool,
            line_no: int,
        ) -> None:
            syns = parse_paren_synonyms(left, current_genus, genera, abbrev_index)
            picked = pick_taxon_raw(left, current_genus, genera, abbrev_index)
            if not picked:
                return
            taxon_raw, genus_hint = picked

            expanded = expand_taxon(
                taxon_raw, genus_hint or current_genus, genera, abbrev_index
            )
            primary = name_for_canonical(expanded)
            canon = canonical_binomial(primary)
            if not canon:
                cleaned = re.sub(r"\([^)]*\)", " ", primary)
                cleaned = re.sub(r"\s+", " ", cleaned).strip()
                canon = canonical_binomial(cleaned)
            if not canon:
                return
            g = canon.split()[0]
            if g not in genera and g != current_genus and g != genus_hint:
                return

            syn_canons = [
                s["synonym_canonical"] for s in syns if s["synonym_canonical"] != canon
            ]
            syn_raws = [s["synonym_raw"] for s in syns if s["synonym_canonical"] != canon]
            # In checklist even without island codes when clave says not naturalized.
            in_mll = ("Mll" in islands) or no_naturalitzada

            if canon in by_canon:
                prev = by_canon[canon]
                prev_islands = set(prev["islands"].split(",")) if prev["islands"] else set()
                prev_islands.update(islands)
                prev["islands"] = ",".join(
                    x for x in ("Mll", "Mn", "Ei", "Fo") if x in prev_islands
                )
                prev["in_mallorca"] = bool(prev.get("in_mallorca")) or in_mll
                if no_naturalitzada:
                    prev["no_naturalitzada"] = True
                    note = str(prev.get("extraction_note", "") or "")
                    if "no_naturalitzada" not in note:
                        prev["extraction_note"] = (
                            (note + "; " if note else "") + "no_naturalitzada"
                        )
                prev_syn = [x for x in str(prev.get("paren_synonyms", "")).split(";") if x]
                for sc in syn_canons:
                    if sc not in prev_syn:
                        prev_syn.append(sc)
                prev["paren_synonyms"] = ";".join(prev_syn)
                prev_raw = [
                    x for x in str(prev.get("paren_synonyms_raw", "")).split(" | ") if x
                ]
                for sr in syn_raws:
                    if sr not in prev_raw:
                        prev_raw.append(sr)
                prev["paren_synonyms_raw"] = " | ".join(prev_raw)
                return

            taxon_clave_out = expanded
            if syn_raws and "(=" not in expanded:
                taxon_clave_out = f"{expanded} (= {'; '.join(syn_raws)})"
            if no_naturalitzada and "(no naturalitzada)" not in taxon_clave_out.lower():
                taxon_clave_out = f"{taxon_clave_out} (no naturalitzada)"

            by_canon[canon] = {
                "taxon_clave_raw": taxon_raw,
                "taxon_clave": taxon_clave_out,
                "canonical": canon,
                "species_key": species_level_key(canon),
                "paren_synonyms": ";".join(syn_canons),
                "paren_synonyms_raw": " | ".join(syn_raws),
                "rank_inferred": (
                    "subspecies"
                    if "subsp." in canon
                    else "variety"
                    if " var." in canon
                    else "form"
                    if " f." in canon
                    else "hybrid"
                    if " x " in canon
                    else "species"
                ),
                "islands": ",".join(x for x in ("Mll", "Mn", "Ei", "Fo") if x in islands),
                "in_mallorca": in_mll,
                "no_naturalitzada": bool(no_naturalitzada),
                "genus_context": genus_hint or current_genus,
                "family_context": current_family,
                "source_line": line_no,
                "extraction_note": "no_naturalitzada" if no_naturalitzada else "",
            }

        for m in ISLAND_RE.finditer(window):
            islands = parse_islands(m.group("islands"))
            if not any(x in islands for x in ("Mll", "Mn", "Ei", "Fo")):
                continue
            if mallorca_only and "Mll" not in islands:
                continue
            left = window[: m.start()].rstrip(" .")
            # Island hit may also carry the non-naturalized note nearby.
            no_nat = bool(NO_NAT_RE.search(window[m.start() : m.end() + 40]))
            _record_from_left(left, islands, no_naturalitzada=no_nat, line_no=i + 1)

        for m in NO_NAT_RE.finditer(window):
            # Skip if this marker sits inside an island-code paren already handled.
            # Island codes never use the no-naturalitzada wording.
            left = window[: m.start()].rstrip(" .")
            if mallorca_only:
                # Still include: clave lists them for the flora even if not naturalized.
                pass
            _record_from_left(left, [], no_naturalitzada=True, line_no=i + 1)

    # Fix mis-parses: "... ebusitanum ... (inclou L. pseudebusitanum)"
    # If we stored the synonym, also ensure primary is present.
    for canon, row in list(by_canon.items()):
        raw = str(row.get("taxon_clave_raw", "")) + " " + str(row.get("taxon_clave", ""))
        if re.search(r"ebusitanum.*inclou|inclou.*pseudebusitanum", raw, re.I) or (
            "pseudebusitanum" in canon and "Limonium" in canon
        ):
            # check source window around line for ebusitanum primary
            pass
    # Prefer Limonium ebusitanum over pseudebusitanum when both from same inclou construct
    if "Limonium pseudebusitanum" in by_canon and "Limonium ebusitanum" not in by_canon:
        row = by_canon.pop("Limonium pseudebusitanum")
        row = dict(row)
        row["canonical"] = "Limonium ebusitanum"
        row["species_key"] = "Limonium ebusitanum"
        row["taxon_clave"] = "Limonium ebusitanum (Font Quer) Font Quer (inclou L. pseudebusitanum)"
        row["extraction_note"] = "renamed_from_inclou_synonym_pseudebusitanum"
        by_canon["Limonium ebusitanum"] = row

    extras = extract_unlabeled_extras(text, genera, abbrev_index, by_canon)
    if extras:
        print(f"added unlabeled extras: {len(extras)}")

    df = pd.DataFrame(list(by_canon.values()))
    if df.empty:
        return df
    if "extraction_note" not in df.columns:
        df["extraction_note"] = ""
    if "paren_synonyms" not in df.columns:
        df["paren_synonyms"] = ""
    if "paren_synonyms_raw" not in df.columns:
        df["paren_synonyms_raw"] = ""
    if "no_naturalitzada" not in df.columns:
        df["no_naturalitzada"] = False
    df["paren_synonyms"] = df["paren_synonyms"].fillna("")
    df["paren_synonyms_raw"] = df["paren_synonyms_raw"].fillna("")
    df["no_naturalitzada"] = df["no_naturalitzada"].fillna(False).astype(bool)
    df = df.sort_values("canonical").reset_index(drop=True)
    df.insert(0, "clave_id", range(1, len(df) + 1))
    return df


def build_dataset_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    accepted = pd.read_excel(DATASET_XLSX, sheet_name="gpkg_accepted")
    variants = pd.read_excel(DATASET_XLSX, sheet_name="gpkg_original_variants")

    accepted = accepted.copy()
    accepted["canonical"] = accepted["accepted_species_name"].map(canonical_binomial)
    accepted["species_key_name"] = accepted["canonical"].map(species_level_key)
    accepted["accepted_species_key"] = pd.to_numeric(
        accepted["accepted_species_key"], errors="coerce"
    ).astype("Int64")

    variants = variants.copy()
    variants["canonical_original"] = variants["original_name"].map(canonical_binomial)
    variants["canonical_accepted"] = variants["accepted_species_name"].map(canonical_binomial)
    return accepted, variants


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


def resolve_taxon_gbif(name: str) -> dict:
    empty = {
        "query_name": name,
        "accepted_key": pd.NA,
        "accepted_rank": "OTHER",
        "accepted_name": pd.NA,
        "accepted_species_key": pd.NA,
        "accepted_species_name": pd.NA,
        "accepted_order": pd.NA,
        "accepted_family": pd.NA,
        "gbif_match_type": pd.NA,
        "gbif_confidence": pd.NA,
        "gbif_status": "not_found",
    }
    match = gbif_get(
        "https://api.gbif.org/v1/species/match?"
        f"name={quote(name)}&kingdom=Plantae"
    )
    if not match or not match.get("usageKey"):
        return empty
    rank = match.get("rank") or "OTHER"
    if rank in {"GENUS", "FAMILY", "CLASS", "ORDER", "PHYLUM", "KINGDOM"}:
        empty["accepted_rank"] = rank
        empty["accepted_order"] = match.get("order")
        empty["accepted_family"] = match.get("family")
        empty["gbif_match_type"] = match.get("matchType")
        empty["gbif_confidence"] = match.get("confidence")
        empty["gbif_status"] = "high_rank_only"
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
        "query_name": name,
        "accepted_key": accepted_key,
        "accepted_rank": rank,
        "accepted_name": accepted_name,
        "accepted_species_key": species_key,
        "accepted_species_name": species_name,
        "accepted_order": match.get("order"),
        "accepted_family": match.get("family"),
        "gbif_match_type": match.get("matchType"),
        "gbif_confidence": match.get("confidence"),
        "gbif_status": "ok" if pd.notna(species_key) else "no_species_key",
    }


def _norm_key(s: str | None) -> str | None:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return None
    s = re.sub(r"\s+", " ", str(s)).strip()
    return s.casefold() if s else None


def load_taxonkey_caches() -> dict[str, dict]:
    """Map normalized query string -> GBIF row dict."""
    by_name: dict[str, dict] = {}
    frames = []
    for path in (TAXONKEY_REVIEW_XLSX, TAXONKEY_XLSX, GBIF_CACHE_XLSX):
        if path.exists():
            df = pd.read_excel(path)
            frames.append(df)
            print(f"cache load {path.name}: {len(df)} rows")
    if not frames:
        return by_name

    for df in frames:
        if "initial_name" in df.columns and "query_name" not in df.columns:
            df = df.rename(columns={"initial_name": "query_name"})
        for _, row in df.iterrows():
            q = row.get("query_name")
            nk = _norm_key(q)
            if not nk:
                continue
            # prefer rows that already have species key
            prev = by_name.get(nk)
            if prev is None or (
                pd.isna(prev.get("accepted_species_key"))
                and pd.notna(row.get("accepted_species_key"))
            ):
                by_name[nk] = row.to_dict()
            # also index by canonical of query / accepted name
            for field in ("query_name", "accepted_name", "accepted_species_name"):
                c = canonical_binomial(str(row.get(field))) if pd.notna(row.get(field)) else None
                if c:
                    ck = _norm_key(c)
                    if ck and ck not in by_name:
                        by_name[ck] = row.to_dict()
                    sk = species_level_key(c)
                    if sk:
                        skk = _norm_key(sk)
                        if skk and skk not in by_name:
                            by_name[skk] = row.to_dict()
    return by_name


def query_names_for_clave_row(row: pd.Series) -> list[str]:
    """Ordered GBIF query candidates for one clave taxon."""
    names: list[str] = []
    primary = name_for_canonical(str(row["taxon_clave"]))
    if primary:
        names.append(primary)
    # drop parenthetical authors for a cleaner binomial+infra query
    no_paren = re.sub(r"\([^)]*\)", " ", primary)
    no_paren = re.sub(r"\s+", " ", no_paren).strip(" .")
    if no_paren and no_paren not in names:
        names.append(no_paren)
    canon = row.get("canonical")
    if pd.notna(canon) and canon not in names:
        names.append(str(canon))
    sp = row.get("species_key")
    if pd.notna(sp) and str(sp) not in names:
        names.append(str(sp))
    # unique preserve order
    out, seen = [], set()
    for n in names:
        k = _norm_key(n)
        if k and k not in seen:
            seen.add(k)
            out.append(n)
    return out


def resolve_clave_gbif(clave: pd.DataFrame, fetch: bool = True) -> pd.DataFrame:
    cache = load_taxonkey_caches()
    rows = []
    fetched = 0
    reused = 0
    new_cache_rows = []

    for i, crow in clave.iterrows():
        queries = query_names_for_clave_row(crow)
        hit = None
        used_query = None
        source = None
        for q in queries:
            nk = _norm_key(q)
            if nk and nk in cache and pd.notna(cache[nk].get("accepted_species_key")):
                hit = cache[nk]
                used_query = q
                source = "cache"
                reused += 1
                break
        if hit is None and fetch:
            for q in queries:
                if (i + 1) % 50 == 0 or fetched == 0:
                    print(f"  GBIF fetch ~{fetched}+ for clave_id={crow['clave_id']} ({q})")
                resolved = resolve_taxon_gbif(q)
                fetched += 1
                time.sleep(0.05)
                # store all attempts in cache map
                cache[_norm_key(q)] = {
                    **resolved,
                    "query_name": q,
                }
                new_cache_rows.append(cache[_norm_key(q)])
                if resolved.get("gbif_status") == "ok" and pd.notna(
                    resolved.get("accepted_species_key")
                ):
                    hit = resolved
                    used_query = q
                    source = "gbif_api"
                    break
        if hit is None:
            # last cache miss without species key, if any
            for q in queries:
                nk = _norm_key(q)
                if nk and nk in cache:
                    hit = cache[nk]
                    used_query = q
                    source = "cache_weak"
                    break

        def _get(h, k, default=pd.NA):
            if h is None:
                return default
            v = h.get(k, default)
            return default if v is None else v

        rows.append(
            {
                **crow.to_dict(),
                "gbif_query": used_query,
                "gbif_source": source or "none",
                "gbif_accepted_key": _get(hit, "accepted_key"),
                "gbif_accepted_rank": _get(hit, "accepted_rank"),
                "gbif_accepted_name": _get(hit, "accepted_name"),
                "gbif_accepted_species_key": _get(hit, "accepted_species_key"),
                "gbif_accepted_species_name": _get(hit, "accepted_species_name"),
                "gbif_order": _get(hit, "accepted_order"),
                "gbif_family": _get(hit, "accepted_family"),
                "gbif_match_type": _get(hit, "gbif_match_type"),
                "gbif_confidence": _get(hit, "gbif_confidence"),
                "gbif_status": _get(hit, "gbif_status", "none"),
            }
        )

    out = pd.DataFrame(rows)
    out["gbif_accepted_species_key"] = pd.to_numeric(
        out["gbif_accepted_species_key"], errors="coerce"
    ).astype("Int64")

    if new_cache_rows:
        cache_df = pd.DataFrame(new_cache_rows)
        if GBIF_CACHE_XLSX.exists():
            old = pd.read_excel(GBIF_CACHE_XLSX)
            cache_df = pd.concat([old, cache_df], ignore_index=True)
        cache_df = cache_df.drop_duplicates("query_name", keep="last")
        cache_df.to_excel(GBIF_CACHE_XLSX, index=False)
        print(f"wrote GBIF cache {GBIF_CACHE_XLSX} rows={len(cache_df)}")

    print(
        f"GBIF resolve: reused_cache_hits={reused} api_calls={fetched} "
        f"with_species_key={(out['gbif_accepted_species_key'].notna()).sum()}"
    )
    return out


def match_clave_to_dataset(clave_gbif: pd.DataFrame, accepted: pd.DataFrame) -> pd.DataFrame:
    acc_by_key: dict[int, int] = {}
    for idx, row in accepted.iterrows():
        k = row["accepted_species_key"]
        if pd.notna(k) and int(k) not in acc_by_key:
            acc_by_key[int(k)] = idx

    acc_by_canon: dict[str, list[int]] = defaultdict(list)
    acc_by_species: dict[str, list[int]] = defaultdict(list)
    for idx, row in accepted.iterrows():
        if pd.notna(row["canonical"]):
            acc_by_canon[row["canonical"]].append(idx)
            sp = species_level_key(row["canonical"])
            if sp:
                acc_by_species[sp].append(idx)

    by_genus: dict[str, list[str]] = defaultdict(list)
    for c in accepted["canonical"].dropna().unique():
        by_genus[c.split()[0].lower()].append(c)

    rows = []
    for _, crow in clave_gbif.iterrows():
        canon = crow["canonical"]
        gbif_key = crow["gbif_accepted_species_key"]
        match_type = "none"
        score = 0.0
        matched_name = None
        matched_canon = None
        matched_key = None
        note = ""
        name_conflict = False

        nearest_name = None
        nearest_canon = None
        nearest_score = None
        nearest_note = ""

        # Name identity (exact or species-level collapse of subsp/var) preferred over
        # GBIF synonym hops that orphan the dataset row (Allium roseum; Ch. album
        # var. album -> berlandieri, etc.).
        name_hit_idx = acc_by_canon[canon][0] if pd.notna(canon) and canon in acc_by_canon else None
        sp_canon = species_level_key(canon) if pd.notna(canon) else None
        sp_hit_idx = None
        if name_hit_idx is None and sp_canon and sp_canon in acc_by_species:
            # Prefer a dataset row whose canonical is exactly the species binomial.
            sp_idxs = acc_by_species[sp_canon]
            exact_sp = [i for i in sp_idxs if accepted.at[i, "canonical"] == sp_canon]
            sp_hit_idx = (exact_sp or sp_idxs)[0]
        gbif_hit_idx = (
            acc_by_key[int(gbif_key)]
            if pd.notna(gbif_key) and int(gbif_key) in acc_by_key
            else None
        )

        name_or_sp_idx = name_hit_idx if name_hit_idx is not None else sp_hit_idx

        # If name/species-level and GBIF disagree:
        # - Exact canonical match always wins (Arundo donax L. vs Arundo donax L;
        #   GBIF wrongly maps donax -> Phragmites).
        # - Species-level collapse only (subsp/var): trust GBIF when genus changes
        #   (Ranunculus ficaria subsp. -> Ficaria verna); trust name when same genus
        #   / other epithet (Ch. album var. album -> berlandieri).
        prefer_gbif_over_name = False
        if name_or_sp_idx is not None and gbif_hit_idx is not None:
            name_ds_key = accepted.at[name_or_sp_idx, "accepted_species_key"]
            gbif_ds_key = accepted.at[gbif_hit_idx, "accepted_species_key"]
            if pd.notna(name_ds_key) and pd.notna(gbif_ds_key) and int(name_ds_key) != int(gbif_ds_key):
                if name_hit_idx is None:
                    gbif_acc = crow.get("gbif_accepted_species_name") or crow.get("gbif_accepted_name")
                    gbif_acc_canon = canonical_binomial(str(gbif_acc)) if pd.notna(gbif_acc) else None
                    clave_genus = str(sp_canon or canon).split()[0].lower() if (sp_canon or canon) else ""
                    gbif_genus = gbif_acc_canon.split()[0].lower() if gbif_acc_canon else ""
                    if clave_genus and gbif_genus and clave_genus != gbif_genus:
                        prefer_gbif_over_name = True

        if name_or_sp_idx is not None and not prefer_gbif_over_name:
            matched_name = accepted.at[name_or_sp_idx, "accepted_species_name"]
            matched_canon = accepted.at[name_or_sp_idx, "canonical"]
            matched_key = accepted.at[name_or_sp_idx, "accepted_species_key"]
            score = 100.0
            via_sp = name_hit_idx is None
            if gbif_hit_idx is not None:
                gbif_ds_key = accepted.at[gbif_hit_idx, "accepted_species_key"]
                if pd.notna(gbif_ds_key) and int(gbif_ds_key) == int(matched_key):
                    match_type = "gbif_species_key"
                    note = f"same accepted_species_key={int(gbif_key)}"
                else:
                    match_type = "name_exact_species_level" if via_sp else "name_exact"
                    gbif_canon = accepted.at[gbif_hit_idx, "canonical"]
                    if via_sp:
                        note = (
                            f"species-level identity '{sp_canon}' "
                            f"(clave '{canon}') preferred over GBIF synonym hop "
                            f"to '{gbif_canon}' (key={int(gbif_key)})"
                        )
                    else:
                        note = (
                            f"name identity '{canon}' preferred over GBIF synonym hop "
                            f"to '{gbif_canon}' (key={int(gbif_key)})"
                        )
            elif via_sp:
                match_type = "name_exact_species_level"
                note = (
                    f"species-level identity '{sp_canon}' from clave '{canon}'"
                    + (
                        "; no GBIF species key"
                        if pd.isna(gbif_key)
                        else f"; GBIF key={int(gbif_key)} not in dataset"
                    )
                )
            elif pd.isna(gbif_key):
                match_type = "name_exact"
                note = "name identity; no GBIF species key"
            else:
                match_type = "name_exact"
                note = f"name identity; GBIF key={int(gbif_key)} not in dataset"
        elif gbif_hit_idx is not None:
            match_type = "gbif_species_key"
            score = 100.0
            matched_name = accepted.at[gbif_hit_idx, "accepted_species_name"]
            matched_canon = accepted.at[gbif_hit_idx, "canonical"]
            matched_key = accepted.at[gbif_hit_idx, "accepted_species_key"]
            if prefer_gbif_over_name and name_or_sp_idx is not None:
                name_canon = accepted.at[name_or_sp_idx, "canonical"]
                note = (
                    f"GBIF genus-change synonymy: clave '{canon}' == dataset "
                    f"'{matched_canon}' (key={int(gbif_key)}); preferred over "
                    f"homonym/species-level name '{name_canon}'"
                )
            elif (
                matched_canon
                and canon
                and matched_canon.split()[0].lower() != str(canon).split()[0].lower()
            ):
                note = (
                    f"GBIF synonymy: clave '{canon}' == dataset '{matched_canon}' "
                    f"(key={int(gbif_key)})"
                )
            else:
                note = f"same accepted_species_key={int(gbif_key)}"
        else:
            genus = str(canon).split()[0].lower() if pd.notna(canon) else ""
            pool = by_genus.get(genus, [])
            best = process.extractOne(str(canon), pool, scorer=fuzz.token_sort_ratio) if pool else None
            if best and best[1] >= 90:
                matched_canon, score, _ = best[0], float(best[1]), best[2]
                idx = acc_by_canon[matched_canon][0]
                matched_name = accepted.at[idx, "accepted_species_name"]
                matched_key = accepted.at[idx, "accepted_species_key"]
                match_type = "fuzzy_name_fallback"
                note = "name fuzzy fallback; review (no GBIF key join)"
            elif best:
                nearest_canon, nearest_score, _ = best[0], float(best[1]), best[2]
                idx = acc_by_canon[nearest_canon][0]
                nearest_name = accepted.at[idx, "accepted_species_name"]
                nearest_note = "nearest same-genus name (not a match)"
                note = (
                    f"GBIF key={gbif_key} not in dataset"
                    if pd.notna(gbif_key)
                    else "no GBIF key and no name match"
                )
            else:
                note = (
                    f"GBIF key={gbif_key} not in dataset"
                    if pd.notna(gbif_key)
                    else "no GBIF key and no name match"
                )

        # Only fuzzy name joins need manual review. Exact binomial identity does not.
        dubious = match_type == "fuzzy_name_fallback"
        if pd.notna(crow.get("gbif_confidence")) and crow["gbif_confidence"] < 90 and match_type == "gbif_species_key":
            # low confidence GBIF backbone hit
            if crow.get("gbif_match_type") not in {"EXACT", "EXACT_CANONICAL_NAME"}:
                dubious = True
                note = (note + "; " if note else "") + "low/fuzzy GBIF match confidence"

        rows.append(
            {
                "clave_id": crow["clave_id"],
                "taxon_clave": crow["taxon_clave"],
                "canonical_clave": canon,
                "species_key_clave": crow["species_key"],
                "rank_clave": crow["rank_inferred"],
                "islands": crow["islands"],
                "genus_context": crow["genus_context"],
                "gbif_query": crow["gbif_query"],
                "gbif_source": crow["gbif_source"],
                "gbif_accepted_species_key": gbif_key,
                "gbif_accepted_species_name": crow["gbif_accepted_species_name"],
                "gbif_accepted_name": crow["gbif_accepted_name"],
                "gbif_match_type": crow["gbif_match_type"],
                "gbif_confidence": crow["gbif_confidence"],
                "gbif_status": crow["gbif_status"],
                "match_type": match_type,
                "match_score": score,
                "dubious_hook": dubious,
                "name_gbif_conflict": name_conflict,
                "match_note": note,
                "no_naturalitzada": bool(crow.get("no_naturalitzada", False)),
                "dataset_accepted_name": matched_name,
                "dataset_canonical": matched_canon,
                "dataset_accepted_species_key": matched_key,
                "nearest_dataset_name": nearest_name,
                "nearest_dataset_canonical": nearest_canon,
                "nearest_score": nearest_score,
                "nearest_note": nearest_note,
                "exclude_from_dataset_manual": "",
                "review_comment": "",
            }
        )

        # Also hook dataset rows that match (= synonym) names in the clave text.
        primary_ds_key = int(matched_key) if pd.notna(matched_key) else None
        paren = str(crow.get("paren_synonyms", "") or "")
        for syn_canon in [x.strip() for x in paren.split(";") if x.strip()]:
            syn_idx = None
            if syn_canon in acc_by_canon:
                syn_idx = acc_by_canon[syn_canon][0]
            else:
                syn_sp = species_level_key(syn_canon)
                if syn_sp and syn_sp in acc_by_species:
                    exact_sp = [
                        i
                        for i in acc_by_species[syn_sp]
                        if accepted.at[i, "canonical"] == syn_sp
                    ]
                    syn_idx = (exact_sp or acc_by_species[syn_sp])[0]
            if syn_idx is None:
                continue
            syn_ds_key = accepted.at[syn_idx, "accepted_species_key"]
            if primary_ds_key is not None and pd.notna(syn_ds_key) and int(syn_ds_key) == primary_ds_key:
                continue
            rows.append(
                {
                    "clave_id": crow["clave_id"],
                    "taxon_clave": crow["taxon_clave"],
                    "canonical_clave": canon,
                    "species_key_clave": crow["species_key"],
                    "rank_clave": crow["rank_inferred"],
                    "islands": crow["islands"],
                    "genus_context": crow["genus_context"],
                    "gbif_query": crow["gbif_query"],
                    "gbif_source": crow["gbif_source"],
                    "gbif_accepted_species_key": gbif_key,
                    "gbif_accepted_species_name": crow["gbif_accepted_species_name"],
                    "gbif_accepted_name": crow["gbif_accepted_name"],
                    "gbif_match_type": crow["gbif_match_type"],
                    "gbif_confidence": crow["gbif_confidence"],
                    "gbif_status": crow["gbif_status"],
                    "match_type": "clave_paren_synonym",
                    "match_score": 100.0,
                    "dubious_hook": False,
                    "name_gbif_conflict": False,
                    "match_note": (
                        f"dataset matches clave paren synonym '{syn_canon}' "
                        f"of primary '{canon}'"
                    ),
                    "no_naturalitzada": bool(crow.get("no_naturalitzada", False)),
                    "dataset_accepted_name": accepted.at[syn_idx, "accepted_species_name"],
                    "dataset_canonical": accepted.at[syn_idx, "canonical"],
                    "dataset_accepted_species_key": syn_ds_key,
                    "nearest_dataset_name": None,
                    "nearest_dataset_canonical": None,
                    "nearest_score": None,
                    "nearest_note": "",
                    "exclude_from_dataset_manual": "",
                    "review_comment": "",
                }
            )

    return pd.DataFrame(rows)


def main() -> None:
    text, source_label = load_clave_text()
    print(f"Reading clave source: {source_label}")
    clave_all = extract_clave_taxa(text, mallorca_only=False)
    clave = clave_all[clave_all["in_mallorca"]].copy()
    clave = clave.reset_index(drop=True)
    clave["clave_id"] = range(1, len(clave) + 1)
    print(
        f"Clave taxa extracted: all_islands={len(clave_all)} "
        f"Mallorca(Mll)={len(clave)} "
        f"(Quarta edició claims 2056; old PDF claimed 1975)"
    )

    accepted, _variants = build_dataset_tables()
    print(f"Dataset accepted species: {len(accepted)}")

    # GBIF resolve for ALL clave taxa (cache shared); match dataset on Mll subset
    clave_all_gbif = resolve_clave_gbif(clave_all, fetch=True)
    clave_gbif = clave_all_gbif[clave_all_gbif["in_mallorca"]].copy()
    clave_gbif = clave_gbif.reset_index(drop=True)
    clave_gbif["clave_id"] = range(1, len(clave_gbif) + 1)
    matches = match_clave_to_dataset(clave_gbif, accepted)

    solid_types = {
        "gbif_species_key",
        "name_exact",
        "name_exact_species_level",
        "name_exact_no_gbif",
        "name_exact_gbif_miss",
        "clave_paren_synonym",
    }
    hooked_types = solid_types | {"fuzzy_name_fallback"}
    matched_ok = matches[matches["match_type"].isin(solid_types)]
    matched_dubious = matches[matches["dubious_hook"] == True]  # noqa: E712
    unmatched_clave = matches[matches["match_type"].eq("none")]

    # Dataset row is "in clave" if a successful match pointed at its key
    # (name identity or GBIF), not merely if clave's GBIF synonym key equals it.
    hooked_mask = matches["match_type"].isin(hooked_types)
    matched_dataset_keys = set(
        int(k)
        for k in matches.loc[hooked_mask, "dataset_accepted_species_key"].dropna().unique()
    )

    accepted2 = accepted.copy()
    accepted2["clave_status"] = accepted2["accepted_species_key"].map(
        lambda k: "present_in_clave"
        if pd.notna(k) and int(k) in matched_dataset_keys
        else "missing_in_clave"
    )
    dataset_not_in_clave = accepted2[accepted2["clave_status"] == "missing_in_clave"][
        [
            "accepted_species_name",
            "canonical",
            "species_key_name",
            "accepted_species_key",
            "n_obs_gpkg",
            "PAC",
            "PAU",
            "POV",
            "POU",
            "clave_status",
        ]
    ].copy()
    dataset_not_in_clave["exclude_from_dataset_manual"] = ""
    dataset_not_in_clave["review_comment"] = ""
    dataset_not_in_clave = dataset_not_in_clave.sort_values("accepted_species_name")

    clave_not_in_dataset = unmatched_clave.copy()
    clave_not_in_dataset["exclude_or_add_manual"] = ""

    synonym_hooks = matches[
        matches["match_type"].eq("gbif_species_key")
        & matches["match_note"].fillna("").str.contains("GBIF synonymy", regex=False)
    ]

    exclude_candidates = dataset_not_in_clave.copy()
    exclude_candidates.insert(0, "reason", "in_dataset_not_in_clave_gbif_key")

    summary = pd.DataFrame(
        [
            {"metric": "source_file", "value": source_label},
            {"metric": "quarta_claimed_taxa", "value": 2056},
            {"metric": "old_pdf_claimed_taxa", "value": 1975},
            {"metric": "clave_all_islands_extracted", "value": len(clave_all_gbif)},
            {"metric": "clave_Mll_extracted", "value": len(clave_gbif)},
            {"metric": "gap_vs_quarta_2056", "value": 2056 - len(clave_all_gbif)},
            {"metric": "clave_Mll_with_gbif_species_key", "value": int(matches["gbif_accepted_species_key"].notna().sum())},
            {"metric": "dataset_accepted_species", "value": len(accepted)},
            {"metric": "matches_ok", "value": len(matched_ok)},
            {"metric": "matches_gbif_species_key", "value": int(matches["match_type"].eq("gbif_species_key").sum())},
            {"metric": "matches_name_exact", "value": int(matches["match_type"].eq("name_exact").sum())},
            {"metric": "matches_name_exact_species_level", "value": int(matches["match_type"].eq("name_exact_species_level").sum())},
            {"metric": "matches_clave_paren_synonym", "value": int(matches["match_type"].eq("clave_paren_synonym").sum())},
            {"metric": "matches_gbif_synonym_genus_change", "value": len(synonym_hooks)},
            {"metric": "matches_fuzzy_name", "value": int(matches["match_type"].eq("fuzzy_name_fallback").sum())},
            {"metric": "clave_not_in_dataset", "value": len(clave_not_in_dataset)},
            {"metric": "dubious_hooks", "value": int(matches["dubious_hook"].sum())},
            {"metric": "name_gbif_conflicts", "value": int(matches["name_gbif_conflict"].sum())},
            {"metric": "dataset_not_in_clave", "value": len(dataset_not_in_clave)},
            {"metric": "source_dataset", "value": str(DATASET_XLSX.name) + "#gpkg_accepted"},
            {"metric": "join_field", "value": "name_identity_then_accepted_species_key"},
        ]
    )

    clave_out = clave_gbif.copy()
    clave_out["exclude_manual"] = ""
    clave_all_out = clave_all_gbif.copy()
    clave_all_out["exclude_manual"] = ""

    with pd.ExcelWriter(CLAVE_ONLY_XLSX, engine="openpyxl") as xl:
        clave_all_out.to_excel(xl, sheet_name="clave_all", index=False)
        clave_out.to_excel(xl, sheet_name="clave_Mll", index=False)
        summary.iloc[:6].to_excel(xl, sheet_name="summary", index=False)

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xl:
        summary.to_excel(xl, sheet_name="summary", index=False)
        clave_all_out.to_excel(xl, sheet_name="clave_all", index=False)
        clave_out.to_excel(xl, sheet_name="clave_Mll", index=False)
        matches.sort_values(
            ["dubious_hook", "name_gbif_conflict", "match_type", "canonical_clave"],
            ascending=[False, False, True, True],
        ).to_excel(xl, sheet_name="matches", index=False)
        matched_dubious.sort_values(["match_type", "match_score"]).to_excel(
            xl, sheet_name="enganches_dudosos", index=False
        )
        synonym_hooks.sort_values("canonical_clave").to_excel(
            xl, sheet_name="sinonimos_gbif", index=False
        )
        dataset_not_in_clave.to_excel(xl, sheet_name="dataset_no_en_clave", index=False)
        clave_not_in_dataset.to_excel(xl, sheet_name="clave_no_en_dataset", index=False)
        exclude_candidates.to_excel(xl, sheet_name="excluir_candidatos", index=False)

    print(f"Wrote {CLAVE_ONLY_XLSX}")
    print(f"Wrote {OUT_XLSX}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
