#!/usr/bin/env python3
"""Step 1: pull the eligible gene list out of one saved R03 announcement page.

Run from the dark_proteins folder:
    python3 scripts/01_extract_foa_genes.py data/raw/foa_pages/RFA-RM-18-021.html

What it does, in order:
  1. Turns the saved web page into plain text, one line per paragraph/heading.
  2. Finds the line "IDG-eligible proteins open for study under this FOA".
  3. After it, finds the three family headings: Kinases, Ion Channels, GPCRs.
  4. Takes the text under each heading (up to the next heading, or up to the
     normal text that follows the GPCR list) and splits it on commas.
  5. Keeps every piece that looks like a gene symbol, remembering its family.
  6. Prints counts per family and anything unusual, so you can check it.

Output: data/processed/foa_genes/foa_genes_<FOA>.tsv with columns
    foa, family, position, symbol
(position = order in which the symbol appears within its family on the page)
Nothing is renamed or corrected here; symbols are kept exactly as printed.
"""
import csv
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

START = "IDG-eligible proteins open for study"
FAMILY_HEADINGS = {           # how a heading line may be written -> family name we store
    "kinases": "Kinase",
    "ion channels": "Ion channel",
    "gpcrs": "GPCR",
}
SYMBOL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.\-]*$")   # e.g. ADCK1, CACNA2D2, ADGRB3, C1orf167


class TextLines(HTMLParser):
    """Collect visible text, starting a new line at each block-level tag."""
    BLOCK = {"p", "div", "br", "li", "ul", "ol", "h1", "h2", "h3", "h4", "h5", "h6", "tr", "table", "section"}

    def __init__(self):
        super().__init__()
        self.lines, self.buf, self.skip = [], [], 0

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style"):
            self.skip += 1
        if tag in self.BLOCK:
            self.flush()

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self.skip -= 1
        if tag in self.BLOCK:
            self.flush()

    def handle_data(self, data):
        if not self.skip:
            self.buf.append(data)

    def flush(self):
        line = " ".join("".join(self.buf).split())
        if line:
            self.lines.append(line)
        self.buf = []


def page_lines(path):
    p = TextLines()
    p.feed(Path(path).read_text(encoding="utf-8", errors="replace"))
    p.flush()
    return p.lines


def heading_family(line):
    """Return the family if this line is just a family heading (e.g. 'Kinases' or 'Kinases:')."""
    key = line.strip().rstrip(":").strip().lower()
    return FAMILY_HEADINGS.get(key)


def looks_like_prose(line):
    """Normal sentences have several lower-case words; gene lists do not."""
    return len(re.findall(r"\b[a-z]{3,}\b", line)) >= 4


def extract(lines):
    try:
        start = next(i for i, l in enumerate(lines) if START.lower() in l.lower())
    except StopIteration:
        raise SystemExit(f"Could not find the line '{START}' in the page.")
    families, current, odd, marked = {}, None, [], []
    for line in lines[start + 1:]:
        fam = heading_family(line)
        if fam:
            current = fam
            families.setdefault(fam, [])
            continue
        if current is None:
            continue
        if looks_like_prose(line):
            if current == "GPCR":
                break               # normal text after the last list: we are done
            continue
        for piece in line.split(","):
            raw = piece.strip()
            token = raw.strip("*").strip()
            if not token:
                continue
            if token != raw:
                marked.append((current, raw))   # e.g. a footnote star on the page
            if SYMBOL.match(token):
                families[current].append(token)
            else:
                odd.append((current, token))
    return start, families, odd, marked


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 scripts/01_extract_foa_genes.py data/raw/foa_pages/RFA-RM-18-021.html")
    path = Path(sys.argv[1])
    foa = path.stem                                   # e.g. RFA-RM-18-021
    lines = page_lines(path)
    start, families, odd, marked = extract(lines)

    print(f"Page: {path.name}  ({len(lines)} text lines)")
    print(f"Found the list heading on line {start + 1}: \"{lines[start][:80]}\"")
    print(f"Families found: {list(families)}")
    total = 0
    for fam, syms in families.items():
        dups = sorted({s for s in syms if syms.count(s) > 1})
        total += len(syms)
        print(f"\n{fam}: {len(syms)} symbols"
              + (f"  | appears twice: {', '.join(dups)}" if dups else ""))
        print(f"   first 5: {', '.join(syms[:5])}")
        print(f"   last 5:  {', '.join(syms[-5:])}")
    across = {}
    for fam, syms in families.items():
        for s in set(syms):
            across.setdefault(s, []).append(fam)
    multi = {s: f for s, f in across.items() if len(f) > 1}
    print(f"\nTotal: {total} symbols, {len(across)} distinct")
    print(f"In more than one family: {multi or 'none'}")
    print(f"Pieces that did not look like a gene symbol: {odd or 'none'}")
    print(f"Symbols printed with a mark like * (mark removed, check the page): {marked or 'none'}")

    out = Path("data/processed/foa_genes") / f"foa_genes_{foa}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["foa", "family", "position", "symbol"])
        for fam, syms in families.items():
            for i, s in enumerate(syms, 1):
                w.writerow([foa, fam, i, s])
    print(f"\nSaved: {out}")


if __name__ == "__main__":
    main()
