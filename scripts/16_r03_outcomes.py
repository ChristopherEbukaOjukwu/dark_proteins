#!/usr/bin/env python3
"""Analysis 6: what the R03 awards reported. One row per award (all 98), from three sources.

Run from the dark_proteins folder, after scripts 03 and 05:
    python3 scripts/16_r03_outcomes.py --email you@example.edu
It asks NIH RePORTER for each award's "Project Outcomes" text (about 100 requests, ~2 min) and NCBI for the
PMIDs and years of the papers IDG lists; all answers are saved and reused, so reruns are offline.

Sources, per award
    1  RePORTER Project Outcomes    the PI's end-of-grant summary, verbatim
                                    (reporter.nih.gov/services/Projects/ProjectOutcomes?projectId=<appl_id>;
                                    every appl_id of the award is asked, the first non-empty text is kept)
    2  IDG award list               the "Outcomes/Publications" column of R03_Award_List.csv: free text and/or
                                    paper links (PMC and PubMed links turned into PMIDs)
    3  RePORTER-linked papers       papers that list the grant number (script 05, r03_papers.tsv), with status
For IDG-listed papers the same tie rule as script 05: about a funded gene = the paper mentions one of the
award's confirmed genes (gene_papers.tsv, any species), and published in or after the award's start year.

Scope (label, nothing dropped)
    2018 list            the award funded at least one gene from the 2018 list (83)
    2019 newcomers only  it funded only genes first listed in 2019 (4)
    no confirmed gene    family-level or unclear (11)

Outputs
    data/processed/r03_awards/r03_award_outcomes.tsv   one row per award: sources, counts, outcome text
    data/raw/r03_awards/r03_outcome_labels.csv          for you: one row per award with the RePORTER outcome text,
                                                        any IDG prose (15+ words once links are removed) and empty
                                                        label columns; existing labels are never overwritten
Saved along the way
    data/raw/r03_awards/r03_outcomes_reporter.jsonl     appl_id -> outcome text
    data/raw/pubmed/pmc_to_pmid.tsv                     PMCID -> PMID
    data/raw/pubmed/pmid_years.tsv                      PMID -> year (shared with script 05)
"""
import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

AWARD_LIST = Path("data/raw/r03_awards/R03_Award_List.csv")
REPORTER_JSON = Path("data/raw/r03_awards/reporter_r03_projects.json")
AWARDS = Path("data/processed/r03_awards/foa_genes_with_awards.tsv")
GENES = Path("data/processed/r03_awards/gene_summary.tsv")
GENE_PAPERS = Path("data/processed/attention/gene_papers.tsv")
R03_PAPERS = Path("data/processed/attention/r03_papers.tsv")
OUT = Path("data/processed/r03_awards/r03_award_outcomes.tsv")
LABELS = Path("data/raw/r03_awards/r03_outcome_labels.csv")
CACHE = Path("data/raw/r03_awards/r03_outcomes_reporter.jsonl")
PMC_CACHE = Path("data/raw/pubmed/pmc_to_pmid.tsv")
YEARS = Path("data/raw/pubmed/pmid_years.tsv")
OUTCOMES_URL = "https://reporter.nih.gov/services/Projects/ProjectOutcomes?projectId={}"
IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
LABEL_COLS = ["result_on_gene", "tool_made", "led_to_more_work", "nothing_or_negative", "notes"]


def read_tsv(path):
    if not path.is_file():
        raise SystemExit(f"Missing {path}")
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def core_of(pn):
    m = re.search(r"[A-Z]{2}\d{6}", pn or "")
    return "R03" + m.group(0) if m else ""


def get_json(url, data=None, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "Mozilla/5.0 (dark_proteins research)",
                                                                  "Accept": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=60))
        except Exception as e:
            print(f"   retry {attempt + 1}: {e}")
            time.sleep(4 * (attempt + 1))
    return None


# ---------------------------------------------------------------- 1. RePORTER outcomes
def reporter_outcomes(appl_ids):
    cache = {}
    if CACHE.is_file():
        for line in CACHE.open(encoding="utf-8"):
            d = json.loads(line)
            cache[str(d["appl_id"])] = d["text"]
    todo = [a for a in sorted(appl_ids) if a not in cache]
    if todo:
        print(f"Asking RePORTER for {len(todo)} outcome texts (about {len(todo) * 1.2 / 60:.0f} min) ...")
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8") as out:
        for k, a in enumerate(todo, 1):
            d = get_json(OUTCOMES_URL.format(a))
            if d is None:
                print(f"   no answer for appl_id {a}; it will be asked again next run")
                continue
            text = (d.get("outcomes_text") or "").strip() if isinstance(d, dict) else ""
            cache[a] = text
            out.write(json.dumps({"appl_id": a, "text": text}) + "\n")
            out.flush()
            if k % 20 == 0:
                print(f"   {k} of {len(todo)}")
            time.sleep(1.1)
    return cache


# ---------------------------------------------------------------- 2. IDG-listed papers -> PMIDs and years
def pmc_to_pmid(pmcids, email):
    known = {}
    if PMC_CACHE.is_file():
        known = {r["pmcid"]: r["pmid"] for r in read_tsv(PMC_CACHE)}
    todo = sorted(set(pmcids) - set(known))
    new_file = not PMC_CACHE.is_file()
    PMC_CACHE.parent.mkdir(parents=True, exist_ok=True)
    with PMC_CACHE.open("a", encoding="utf-8") as out:
        if new_file:
            out.write("pmcid\tpmid\n")
        for b in range(0, len(todo), 150):
            batch = todo[b:b + 150]
            q = urllib.parse.urlencode({"ids": ",".join(batch), "format": "json", "tool": "dark_proteins",
                                        "email": email})
            d = get_json(f"{IDCONV}?{q}")
            if d is None:
                raise SystemExit("NCBI ID converter did not answer; rerun later (finished batches are kept).")
            got = {str(r.get("pmcid", "")): str(r.get("pmid") or "") for r in d.get("records", [])}   # NCBI may send numbers
            for p in batch:
                known[p] = got.get(p, "")
                out.write(f"{p}\t{known[p]}\n")
            time.sleep(0.4)
    return known


def pmid_years(pmids, email):
    known = {}
    if YEARS.is_file():
        known = {r["pmid"]: r["year"] for r in read_tsv(YEARS)}
    todo = sorted({str(p) for p in pmids if p} - set(known), key=int)
    if not todo:
        return known
    new_file = not YEARS.is_file()
    with YEARS.open("a", encoding="utf-8") as out:
        if new_file:
            out.write("pmid\tyear\n")
        for b in range(0, len(todo), 200):
            batch = todo[b:b + 200]
            data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(batch), "retmode": "json",
                                           "tool": "dark_proteins", "email": email}).encode()
            d = get_json(ESUMMARY, data=data)
            if d is None:
                raise SystemExit("PubMed did not answer; rerun later.")
            for p in batch:
                rec = d.get("result", {}).get(p, {})
                m = re.match(r"(\d{4})", rec.get("pubdate", "")) or re.match(r"(\d{4})", rec.get("sortpubdate", ""))
                known[p] = m.group(1) if m else ""
                out.write(f"{p}\t{known[p]}\n")
            time.sleep(0.4)
    return known


def links_in(text):
    pmc = re.findall(r"PMC\d+", text or "", flags=re.I)
    pmids = re.findall(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)|ncbi\.nlm\.nih\.gov/pubmed/(\d+)", text or "")
    return [p.upper() for p in pmc], [a or b for a, b in pmids]


def free_text(text):
    """The IDG entry with links and the paper titles that follow them removed (rough): what is left is prose."""
    t = re.sub(r"https?://\S+", " ", text or "")
    return re.sub(r"\s+", " ", t).strip()


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--email", required=True, help="for NCBI's ID converter and PubMed (their rule)")
    a = ap.parse_args()

    listed_2018 = {g["hgnc_id"] for g in read_tsv(GENES) if "2018" in g["listed_years"].split(";")}
    first_listed = {g["hgnc_id"]: g["listed_years"].split(";")[0] for g in read_tsv(GENES)}
    award_genes, award_round = defaultdict(dict), {}
    for r in read_tsv(AWARDS):
        if r["project_number"]:
            c = core_of(r["project_number"])
            award_genes[c][r["hgnc_id"]] = r["gene"]
            award_round[c] = min(award_round.get(c, "9999"), r["year"])

    rep = json.loads(REPORTER_JSON.read_text(encoding="utf-8"))["results"]
    appl, start = defaultdict(set), {}
    for r in rep:
        c = r.get("core_project_num") or core_of(r.get("project_num", ""))
        if r.get("appl_id"):
            appl[c].add(str(r["appl_id"]))
        sd = (r.get("project_start_date") or "")[:4]
        if sd:
            start[c] = min(start.get(c, "9999"), sd)

    with AWARD_LIST.open(encoding="utf-8-sig", newline="") as f:
        idg = list(csv.DictReader(f))
    col = "Outcomes/Publications"
    awards = {}
    for r in idg:
        c = core_of(r["Project number"])
        rnd = str(2017 + int(r["Round of R03"])) if r["Round of R03"].strip().isdigit() else r["Round of R03"]
        awards[c] = {"round": award_round.get(c, rnd), "project_number": r["Project number"],
                     "pi": r["PI Name (Last name, First name)"], "title": r["R03 Title"],
                     "idg_entry": (r.get(col) or "").strip()}
    print(f"Awards on IDG's list: {len(awards)}")

    # 1. RePORTER outcomes
    texts = reporter_outcomes({x for c in awards for x in appl.get(c, ())})

    # 2. IDG-listed papers
    pmc_all, pm_all = set(), set()
    for c, d in awards.items():
        pmc, pm = links_in(d["idg_entry"])
        d["idg_pmc"], d["idg_pm"] = pmc, pm
        pmc_all.update(pmc)
        pm_all.update(pm)
    conv = pmc_to_pmid(pmc_all, a.email) if pmc_all else {}
    for d in awards.values():
        d["idg_pmids"] = sorted(({str(conv.get(p, "")) for p in d["idg_pmc"]} | set(d["idg_pm"])) - {""}, key=int)
        d["idg_unconverted"] = sum(1 for p in d["idg_pmc"] if not conv.get(p))
    years = pmid_years({p for d in awards.values() for p in d["idg_pmids"]}, a.email)
    want = {p for d in awards.values() for p in d["idg_pmids"]}
    mentions = defaultdict(set)                                # pmid -> genes it mentions (any species)
    with GENE_PAPERS.open(encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["pmid"] in want:
                mentions[r["pmid"]].add(r["hgnc_id"])

    # 3. RePORTER-linked papers
    linked = defaultdict(list)
    for r in read_tsv(R03_PAPERS):
        linked[r["award"]].append(r)

    rows = []
    for c, d in sorted(awards.items(), key=lambda kv: (kv[1]["round"], kv[0])):
        hs = set(award_genes.get(c, {}))
        if not hs:
            scope = "no confirmed gene"
        elif hs & listed_2018:
            scope = "2018 list"
        else:
            scope = "2019 newcomers only" if all(first_listed.get(h) == "2019" for h in hs) else "later genes only"
        out_text = next((texts[x] for x in sorted(appl.get(c, ())) if texts.get(x)), "")
        about = [p for p in d["idg_pmids"] if mentions[p] & hs and (years.get(p, "") or "9999") >= start.get(c, "0000")]
        rep_pmids = {r["pmid"] for r in linked.get(c, [])}
        tied = {r["pmid"] for r in linked.get(c, []) if r["status"].startswith("tied")}
        rows.append({
            "award": c, "round": d["round"], "project_number": d["project_number"], "scope": scope,
            "genes": "; ".join(sorted(award_genes.get(c, {}).values())), "pi": d["pi"], "title": d["title"],
            "reporter_outcome": int(bool(out_text)), "reporter_outcome_words": len(out_text.split()),
            "idg_entry": int(bool(d["idg_entry"])), "idg_free_text_words": len(free_text(d["idg_entry"]).split()),
            "idg_papers": len(d["idg_pmids"]), "idg_papers_unconverted": d["idg_unconverted"],
            "idg_papers_about_gene": len(about),
            "idg_papers_not_reporter_linked": len(set(d["idg_pmids"]) - rep_pmids),
            "idg_papers_about_gene_not_tied": len(set(about) - tied),
            "reporter_papers": len(rep_pmids), "reporter_papers_tied": len(tied),
            "any_outcome": int(bool(out_text) or bool(d["idg_entry"]) or bool(rep_pmids)),
            "idg_pmids": "; ".join(d["idg_pmids"]), "idg_pmids_about_gene": "; ".join(about),
            "reporter_outcome_text": out_text, "idg_entry_text": d["idg_entry"],
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    # labels file: keep anything already filled in
    old = {}
    if LABELS.is_file():
        with LABELS.open(encoding="utf-8-sig", newline="") as f:
            old = {r["award"]: r for r in csv.DictReader(f)}
    with LABELS.open("w", encoding="utf-8", newline="") as f:
        cols = ["award", "round", "scope", "genes", "reporter_outcome_text", "idg_list_text"] + LABEL_COLS
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            idg_text = free_text(r["idg_entry_text"]) if r["idg_free_text_words"] >= 15 else ""
            row = {"award": r["award"], "round": r["round"], "scope": r["scope"], "genes": r["genes"],
                   "reporter_outcome_text": r["reporter_outcome_text"], "idg_list_text": idg_text}
            for k in LABEL_COLS:
                row[k] = old.get(r["award"], {}).get(k, "")
            w.writerow(row)

    # ---- summary
    print(f"\nAwards: {len(rows)}")
    for sc in ("2018 list", "2019 newcomers only", "later genes only", "no confirmed gene"):
        g = [r for r in rows if r["scope"] == sc]
        if not g:
            continue
        n = len(g)
        s = lambda k: sum(1 for r in g if r[k])
        print(f"\n{sc} ({n} awards)")
        print(f"   RePORTER outcome text:      {s('reporter_outcome'):3d} ({100 * s('reporter_outcome') / n:.0f}%)")
        print(f"   IDG list entry:             {s('idg_entry'):3d} ({100 * s('idg_entry') / n:.0f}%)")
        print(f"   RePORTER-linked papers:     {s('reporter_papers'):3d} ({100 * s('reporter_papers') / n:.0f}%)")
        print(f"   any of the three:           {s('any_outcome'):3d} ({100 * s('any_outcome') / n:.0f}%)")
        ip = sum(r["idg_papers"] for r in g)
        print(f"   IDG-listed papers: {ip}; about a funded gene: {sum(r['idg_papers_about_gene'] for r in g)}; "
              f"not in RePORTER's links: {sum(r['idg_papers_not_reporter_linked'] for r in g)}; "
              f"about a funded gene but not yet tied: {sum(r['idg_papers_about_gene_not_tied'] for r in g)}; "
              f"links not converted: {sum(r['idg_papers_unconverted'] for r in g)}")
        print("   by round (RePORTER outcome / IDG entry / any):")
        for rd in sorted({r["round"] for r in g}):
            h = [r for r in g if r["round"] == rd]
            print(f"      {rd}: {sum(r['reporter_outcome'] for r in h)} / {sum(r['idg_entry'] for r in h)} / "
                  f"{sum(r['any_outcome'] for r in h)} of {len(h)}")
    nl = sum(1 for r in rows if r["reporter_outcome"] or r["idg_free_text_words"] >= 15)
    print(f"\nSaved: {OUT} ({len(rows)} rows)")
    print(f"Saved: {LABELS}  ({nl} awards with outcome text to label; existing labels kept)")


if __name__ == "__main__":
    main()
