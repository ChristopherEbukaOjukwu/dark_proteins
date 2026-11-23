#!/usr/bin/env python3
"""Step 6: every NIH grant that names one of our genes, per gene per fiscal year (all years to FY2026).

Run from the dark_proteins folder, after step 3:
    python3 scripts/06_gene_grants.py

Inputs
    data/processed/r03_awards/gene_summary.tsv            our 388 genes
    data/processed/r03_awards/foa_genes_with_awards.tsv   which genes each R03 funded
    data/raw/r03_awards/reporter_r03_projects.json        the 98 R03s (for their PIs)
    data/raw/r03_awards/R03_Award_List.csv                the 98 R03 project numbers
    data/raw/hgnc/hgnc_complete_set.txt                   symbol, old symbols, aliases, full name

Search terms for each gene (each searched as an exact phrase in grant titles + abstracts)
    symbol         the current HGNC symbol (always searched)
    full name      the HGNC approved name, e.g. "G protein-coupled receptor 68"
    old symbol / alias   only if safe: 4+ characters, not "ABC-12"-style, not shared with another
                         gene, not another gene's current symbol
    protein name   a few names we met while checking awards (MANUAL_TERMS), e.g. "Kir4.2"
    family name    a name for a whole family (FAMILY_TERMS), e.g. "MRCK"; credits every member,
                   flagged family_name_only and counted in its own column
    Dropped by hand: aliases that gave false hits (DROP_TERMS). SLACK/SLICK count only when the
    text also mentions a channel; NRK and LTK only when it names the kinase, since both are
    also cell-line names (CONTEXT_TERMS).
    A term that returns more than TOO_COMMON grants is dropped as too common (printed).
RePORTER's own hits are then re-checked here: the term must appear as a whole word (symbols) or
as the phrase (names) in the title or abstract. "where" = title, or abstract (title has no match).

Saved along the way (reruns only search new terms)
    data/raw/reporter_grants/term_searches.jsonl     one line per searched term: RePORTER's records

Outputs (data/processed/funding/)
    gene_grants.tsv       one row per gene x grant record (one record = one grant in one fiscal year)
    gene_year_grants.tsv  per gene x fiscal year x grant group: active grants and new grants
                          (new = first fiscal year this grant is found for this gene),
                          all matches and title-only matches (family-name-only hits in their own column)
grant_group: IDG R03 | R01 (R01, R37, RF1) | R21 | R03 (other) | other research | career/training | SBIR/STTR | other
"""
import csv
import json
import re
import time
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

GENES = Path("data/processed/r03_awards/gene_summary.tsv")
AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
R03_PROJECTS = Path("data/raw/r03_awards/reporter_r03_projects.json")
AWARD_LIST = Path("data/raw/r03_awards/R03_Award_List.csv")
HGNC = Path("data/raw/hgnc/hgnc_complete_set.txt")
CACHE = Path("data/raw/reporter_grants/term_searches.jsonl")
OUT_DIR = Path("data/processed/funding")
API = "https://api.reporter.nih.gov/v2/projects/search"
TOO_COMMON = 3000
LAST_FY = 2026
FIELDS = ["ApplId", "ProjectNum", "CoreProjectNum", "ActivityCode", "FiscalYear", "ProjectStartDate",
          "ProjectEndDate", "AwardAmount", "PrincipalInvestigators", "Organization", "AgencyIcAdmin",
          "ProjectTitle", "AbstractText"]
# Protein names met while checking the R03 awards (titles/abstracts used these, not the symbol).
MANUAL_TERMS = {
    "PRKCQ": ["PKC theta", "PKCtheta", "PKC-theta", "protein kinase C theta"],
    "PRKACB": ["PKACbeta", "PKA C beta"],
    "KCNJ15": ["Kir4.2"], "KCNJ16": ["Kir5.1"], "KCNA6": ["Kv1.6"],
    "CDC42BPA": ["MRCKalpha", "MRCK alpha", "MRCK α", "MRCKα"],
    "CDC42BPB": ["MRCKbeta", "MRCK beta", "MRCK β", "MRCKβ"],
    "CDC42BPG": ["MRCKgamma", "MRCK gamma", "MRCK γ", "MRCKγ"],
    "PKD2L2": ["polycystin2-like2", "polycystin-2-like 2", "polycystin-2 like 2"],
    "NPY2R": ["Y2 receptor"], "NPY5R": ["Y5 receptor"],
    "MAP3K10": ["MLK2", "mixed lineage kinase 2"],
    "MAPK4": ["ERK4"], "MAPK15": ["ERK7", "ERK8"],
    "CACNA2D2": ["alpha2delta-2"], "CACNA2D3": ["alpha2delta-3"], "CACNA2D4": ["alpha2delta-4"],
}
# Family names: a hit credits every listed member, labelled "family name" (counted separately).
FAMILY_TERMS = {"MRCK": ["CDC42BPA", "CDC42BPB", "CDC42BPG"], "CDC42BP": ["CDC42BPA", "CDC42BPB", "CDC42BPG"],
                "PCTAIRE": ["CDK16", "CDK17", "CDK18"], "PITSLRE": ["CDK11A", "CDK11B"]}
# Aliases checked by hand and found to give false hits (another gene, an ordinary word, a yeast gene).
DROP_TERMS = {"PAFR", "CSC1", "CDK11", "NIM1",                      # round 1 (top-30 check)
              "CALP", "YAK1", "CDC2L", "MCLC", "PGR2", "IEDA", "LPS1", "LPS2", "BREK", "NERCC",
              "TAPK", "TEIF", "NKTL", "SHIK", "PMLP", "MIDORI", "PSK1", "PSK2", "SLOB", "ERG2",
              "GPR80", "GPHR", "NDR2",                                  # round 2 (full skim)
              "PCTAIRE", "PITSLRE"}                                     # searched as family names instead
# Ordinary words that are also channel names: counted only if the same text mentions a channel.
KNA = r"sodium-activated|Na\+-activated|KNa|KCNT|Slo2|Slack channel|Slick channel|Slack current|Slack and Slick|Slick and Slack"
# Symbols that are also cell-line names: NRK = normal rat kidney cells; "Ltk-" = mouse L cells
# (and "ltk" is a misreading of ITK). Counted only if the same text names the kinase.
NRK_CTX = r"Nik[\s-]*related"
LTK_CTX = r"leukocyte (receptor )?tyrosine kinase|(?<![A-Za-z0-9])ALK(?![A-Za-z0-9])"
CONTEXT_TERMS = {"SLACK": KNA, "SLICK": KNA, "NRK": NRK_CTX, "LTK": LTK_CTX}
HYPHEN_STYLE = re.compile(r"^[A-Za-z]{1,4}-\d{1,3}[A-Za-z]?$")


def read_tsv(path, delimiter="\t"):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_tsv(path, rows, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def core_of(project_number):
    return "R03" + re.search(r"[A-Z]{2}\d{6}", project_number).group(0)


def grant_group(code):
    if code in ("R01", "R37", "RF1"):
        return "R01"
    if code == "R21":
        return "R21"
    if code == "R03":
        return "R03 (other)"
    if code in ("R41", "R42", "R43", "R44", "U43", "U44"):
        return "SBIR/STTR"
    if code[:1] in ("F", "K", "T") or code[:2] == "IK" or code in ("R00", "R25", "R36"):
        return "career/training"
    if code[:1] in ("R", "U", "P", "D", "S", "Z") or code in ("I01", "OT2"):
        return "other research"
    return "other"


# ---------------------------------------------------------------- search terms
def build_terms(ours):
    hg = read_tsv(HGNC)
    current = {r["symbol"].upper() for r in hg}
    alias_owners = defaultdict(set)
    for r in hg:
        for col in ("alias_symbol", "prev_symbol"):
            for a in (r.get(col) or "").replace('"', "").split("|"):
                if a.strip():
                    alias_owners[a.strip().upper()].add(r["hgnc_id"])
    info = {r["hgnc_id"]: r for r in hg}
    terms, skipped = defaultdict(list), Counter()
    for h, sym in ours.items():
        r = info.get(h, {})
        terms[h].append((sym, "symbol"))
        if r.get("name"):
            terms[h].append((r["name"], "full name"))
        for col, label in (("prev_symbol", "old symbol"), ("alias_symbol", "alias")):
            for a in (r.get(col) or "").replace('"', "").split("|"):
                a = a.strip()
                if not a or a.upper() == sym.upper():
                    continue
                if a.upper() in DROP_TERMS:
                    skipped["dropped by hand (false hits)"] += 1
                elif len(a) < 4:
                    skipped["shorter than 4 characters"] += 1
                elif HYPHEN_STYLE.match(a):
                    skipped['"ABC-12"-style'] += 1
                elif a.upper() in current:
                    skipped["another gene's current symbol"] += 1
                elif len(alias_owners[a.upper()]) > 1:
                    skipped["shared by more than one gene"] += 1
                else:
                    terms[h].append((a, label))
        for t in MANUAL_TERMS.get(sym, []):
            terms[h].append((t, "protein name"))
        for t, members in FAMILY_TERMS.items():
            if sym in members:
                terms[h].append((t, "family name"))
        seen, uniq = set(), []
        for t, lab in terms[h]:
            if t.lower() not in seen:
                seen.add(t.lower())
                uniq.append((t, lab))
        terms[h] = uniq
    return terms, skipped


# ---------------------------------------------------------------- RePORTER
def load_cache():
    out = {}
    if CACHE.is_file():
        for line in CACHE.open(encoding="utf-8"):
            d = json.loads(line)
            out[d["term"].lower()] = d
    return out


def search(term, cache, state):
    key = term.lower()
    if key in cache:
        return cache[key]
    results, offset, total = [], 0, None
    while True:
        body = {"criteria": {"advanced_text_search": {"operator": "advanced",
                                                      "search_field": "projecttitle,abstracttext",
                                                      "search_text": f'"{term}"'},
                             "fiscal_years": list(range(1985, LAST_FY + 1))},
                "include_fields": FIELDS, "offset": offset, "limit": 500}
        data = None
        for attempt in range(4):
            wait = 1.1 - (time.monotonic() - state["last"])
            if wait > 0:
                time.sleep(wait)
            try:
                req = urllib.request.Request(API, data=json.dumps(body).encode(),
                                             headers={"Content-Type": "application/json"})
                data = json.load(urllib.request.urlopen(req, timeout=120))
                state["last"] = time.monotonic()
                break
            except Exception as e:
                state["last"] = time.monotonic()
                print(f"   retry {attempt + 1} for {term!r}: {e}")
                time.sleep(5 * (attempt + 1))
        if data is None:
            raise SystemExit("RePORTER did not answer after 4 tries; rerun later (finished terms are kept).")
        total = data.get("meta", {}).get("total", 0)
        if total > TOO_COMMON:
            results = []
            break
        results += data.get("results", [])
        offset += 500
        if offset >= total or offset >= 15000:
            break
    keep = []
    for r in results:
        keep.append({"appl_id": r.get("appl_id"), "project_num": r.get("project_num"),
                     "core": r.get("core_project_num"), "activity": r.get("activity_code") or "",
                     "fy": r.get("fiscal_year"), "start": (r.get("project_start_date") or "")[:10],
                     "end": (r.get("project_end_date") or "")[:10], "amount": r.get("award_amount"),
                     "pis": [(p.get("full_name", ""), p.get("profile_id")) for p in r.get("principal_investigators") or []],
                     "org": (r.get("organization") or {}).get("org_name", ""),
                     "agency": (r.get("agency_ic_admin") or {}).get("abbreviation", ""),
                     "title": r.get("project_title") or "", "abstract": r.get("abstract_text") or ""})
    d = {"term": term, "total": total, "too_common": total > TOO_COMMON, "results": keep}
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8") as f:
        f.write(json.dumps(d) + "\n")
    cache[key] = d
    state["n"] += 1
    return d


def found_in(term, kind, text):
    ctx = CONTEXT_TERMS.get(term.upper())
    if ctx and not re.search(ctx, text, re.I):
        return False
    if kind in ("full name", "protein name") and (" " in term or "-" in term):
        norm = lambda s: re.sub(r"[\s\-]+", " ", s.lower())
        # whole phrase only: "G protein-coupled receptor 6" must not match "... receptor 61"
        return re.search(r"(?<![a-z0-9])" + re.escape(norm(term)) + r"(?![a-z0-9])", norm(text)) is not None
    return re.search(r"(?<![A-Za-z0-9])" + re.escape(term) + r"(?![A-Za-z0-9])", text, re.I) is not None


# ---------------------------------------------------------------- main
def main():
    ours = {r["hgnc_id"]: r["gene"] for r in read_tsv(GENES)}
    terms, skipped = build_terms(ours)
    n_terms = sum(len(v) for v in terms.values())
    print(f"Genes: {len(ours)} | search terms: {n_terms} "
          f"({dict(Counter(lab for v in terms.values() for _, lab in v))})")
    print(f"Aliases skipped as unsafe: {dict(skipped)}")

    # R03s: core numbers, PIs, funded genes
    r03_cores = {core_of(r["Project number"]) for r in read_tsv(AWARD_LIST, ",")}
    r03_pis = defaultdict(set)
    for r in json.loads(R03_PROJECTS.read_text(encoding="utf-8"))["results"]:
        for p in r.get("principal_investigators") or []:
            if p.get("profile_id"):
                r03_pis[r.get("core_project_num")].add(p["profile_id"])
    gene_r03_pis = defaultdict(set)              # hgnc -> PIs of R03s that funded this gene
    for r in read_tsv(AWARDS):
        if r["project_number"]:
            gene_r03_pis[r["hgnc_id"]] |= r03_pis.get(core_of(r["project_number"]), set())

    cache, state = load_cache(), {"last": 0.0, "n": 0}
    todo = sum(1 for v in terms.values() for t, _ in v if t.lower() not in cache)
    print(f"Terms already searched: {n_terms - todo} | to search now: {todo} (about {todo * 1.2 / 60:.0f} min or more)")

    rows, too_common, dropped = {}, [], 0
    for i, (h, tl) in enumerate(sorted(terms.items(), key=lambda kv: ours[kv[0]]), 1):
        for term, kind in tl:
            d = search(term, cache, state)
            if d["too_common"]:
                too_common.append(f"{ours[h]}: {term!r} ({d['total']:,} grants)")
                continue
            for g in d["results"]:
                in_title = found_in(term, kind, g["title"])
                in_abs = found_in(term, kind, g["abstract"])
                if not (in_title or in_abs):
                    dropped += 1
                    continue
                key = (h, g["appl_id"])
                row = rows.get(key)
                if row is None:
                    pis = g["pis"]
                    group = "IDG R03" if g["core"] in r03_cores else grant_group(g["activity"])
                    row = rows[key] = {
                        "gene": ours[h], "hgnc_id": h, "terms": set(), "term_types": set(), "where": "abstract",
                        "appl_id": g["appl_id"], "project_num": g["project_num"], "core_project_num": g["core"],
                        "activity_code": g["activity"], "grant_group": group, "fiscal_year": g["fy"],
                        "project_start": g["start"], "project_end": g["end"], "award_amount": g["amount"],
                        "pi_names": "; ".join(n for n, _ in pis), "pi_profile_ids": "; ".join(str(p) for _, p in pis),
                        "same_pi_as_r03_for_gene": int(bool({p for _, p in pis} & gene_r03_pis.get(h, set()))),
                        "organization": g["org"], "agency": g["agency"], "title": g["title"]}
                row["terms"].add(term)
                row["term_types"].add(kind)
                if in_title:
                    row["where"] = "title"
        if i % 25 == 0:
            print(f"   {i} of {len(terms)} genes done ({state['n']} new searches)")

    out = []
    for r in rows.values():
        out.append({**r, "terms": "; ".join(sorted(r["terms"])), "term_types": "; ".join(sorted(r["term_types"])),
                    "family_name_only": int(r["term_types"] == {"family name"})})
    out.sort(key=lambda r: (r["gene"], r["fiscal_year"] or 0, r["core_project_num"] or ""))

    # per gene x year x group: active and new grants
    first_fy = {}
    for r in out:
        if r["family_name_only"]:          # "new" is judged from gene-specific matches only
            continue
        k = (r["hgnc_id"], r["core_project_num"])
        first_fy[k] = min(first_fy.get(k, 9999), r["fiscal_year"] or 9999)
    agg = defaultdict(lambda: [set(), set(), set(), set(), set()])
    for r in out:
        k = (r["hgnc_id"], r["gene"], r["fiscal_year"], r["grant_group"])
        core = r["core_project_num"]
        if r["family_name_only"]:
            agg[k][4].add(core)
            continue
        new = first_fy[(r["hgnc_id"], core)] == r["fiscal_year"]
        agg[k][0].add(core)
        if new:
            agg[k][1].add(core)
        if r["where"] == "title":
            agg[k][2].add(core)
            if new:
                agg[k][3].add(core)
    yrows = [{"gene": k[1], "hgnc_id": k[0], "fiscal_year": k[2], "grant_group": k[3],
              "n_active": len(v[0]), "n_new": len(v[1]), "n_active_title": len(v[2]), "n_new_title": len(v[3]),
              "n_active_family_name_only": len(v[4])}
             for k, v in sorted(agg.items(), key=lambda kv: (kv[0][1], kv[0][2] or 0, kv[0][3]))]

    # report
    print(f"\nRePORTER hits dropped because the term was not found as a whole word/phrase: {dropped:,}")
    print(f"Terms dropped as too common (> {TOO_COMMON:,} grants): {len(too_common)}")
    for t in too_common:
        print("   " + t)
    print(f"\nGene x grant-year rows: {len(out):,} | distinct grants: {len({r['core_project_num'] for r in out}):,} | "
          f"genes with any grant: {len({r['hgnc_id'] for r in out})} of {len(ours)}")
    print("By grant group (distinct grants, gene-grant pairs):")
    gg = defaultdict(set)
    for r in out:
        gg[r["grant_group"]].add((r["hgnc_id"], r["core_project_num"]))
    for g, s in sorted(gg.items(), key=lambda kv: -len(kv[1])):
        print(f"   {g:<17} {len({c for _, c in s}):>6} grants  {len(s):>6} gene-grant pairs")
    r03_found = {r["core_project_num"] for r in out if r["grant_group"] == "IDG R03"}
    print(f"IDG R03s found by this search: {len(r03_found)} of {len(r03_cores)} (check: should be most)")
    print("\nGrant records per fiscal year, all our genes (to check where coverage starts):")
    fy = Counter(r["fiscal_year"] for r in out)
    for start in range(1985, LAST_FY + 1, 10):
        print("   " + "  ".join(f"{y}:{fy.get(y, 0)}" for y in range(start, min(start + 10, LAST_FY + 1))))
    print("\nCheck three genes (fiscal year: grant groups):")
    for name in ("GPR68", "PRKCQ", "ADGRB3"):
        g = [r for r in out if r["gene"] == name]
        by = defaultdict(set)
        for r in g:
            by[r["fiscal_year"]].add(f"{r['grant_group']}:{r['core_project_num']}")
        print(f"   {name}: " + " | ".join(f"{y}: {', '.join(sorted(v))}" for y, v in sorted(by.items(), key=lambda kv: kv[0] or 0)[-6:]))

    cols = ["gene", "hgnc_id", "terms", "term_types", "family_name_only", "where", "appl_id", "project_num", "core_project_num",
            "activity_code", "grant_group", "fiscal_year", "project_start", "project_end", "award_amount",
            "pi_names", "pi_profile_ids", "same_pi_as_r03_for_gene", "organization", "agency", "title"]
    write_tsv(OUT_DIR / "gene_grants.tsv", out, cols)
    write_tsv(OUT_DIR / "gene_year_grants.tsv", yrows, list(yrows[0]) if yrows else ["gene"])
    print(f"\nSaved: {OUT_DIR / 'gene_grants.tsv'} ({len(out):,} rows)")
    print(f"Saved: {OUT_DIR / 'gene_year_grants.tsv'} ({len(yrows):,} rows)")


if __name__ == "__main__":
    main()
