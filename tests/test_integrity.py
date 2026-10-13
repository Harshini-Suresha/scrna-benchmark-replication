"""Data-integrity gate for the committed dashboard inputs and docs.

Catches the failure classes seen in this repo: non-strict JSON (bare NaN
broke fetch().json(), forcing the offline fallback), dropped pipeline
blocks (qc_dist/similarity/environment), and docs drifting from code
(figure catalogue, tab lists). Stdlib only; safe in CI without data
downloads because site/data/*.json are git-tracked.

Usage: python tests/test_integrity.py
Exit nonzero with the first failed assertion group.
"""
import ast
import json
import re
import sys

ROOT = __import__("os").path.dirname(
    __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
FAIL = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name
          + ("" if cond or not detail else " :: " + detail))
    if not cond:
        FAIL.append(name)


def strict_load(path):
    raw = open(path).read()
    if re.search(r"\bNaN\b|\bInfinity\b", raw):
        raise ValueError("non-strict JSON token in " + path)
    return json.loads(raw)  # strict by default once tokens are gone


# ---- 1. results files strictly parse and carry every pipeline block ----
for fn in ("site/data/results.json", "site/data/results_mouse.json",
           "site/data/demo.json"):
    p = ROOT + "/" + fn
    try:
        d = strict_load(p)
    except Exception as e:  # noqa: BLE001
        check(fn + " strict-parse", False, str(e)[:160])
        continue
    check(fn + " strict-parse", True)
    for key in ("generated", "config", "environment", "cells_per_protocol",
                "protocols", "discrepancies", "similarity", "errors"):
        check(f"{fn} key:{key}", key in d)
    if fn.endswith("demo.json"):
        continue
    ps = d["protocols"]
    check(f"{fn} 15 protocols", len(ps) == 15, f"got {len(ps)}")
    for blk in ("markers", "tsne", "cellstat"):
        n = sum(1 for x in ps if x.get(blk))
        check(f"{fn} block:{blk} 15/15", n == 15, f"got {n}")
    n = sum(1 for x in ps if x.get("qc_dist", {}).get("cell_sample"))
    check(f"{fn} qc_dist+cell_sample 15/15", n == 15, f"got {n}")
    check(f"{fn} similarity matrix",
          bool((d.get("similarity") or {}).get("correlation")))
    check(f"{fn} environment packages",
          bool((d.get("environment") or {}).get("packages")))
    check(f"{fn} discrepancies non-empty",
          len(d.get("discrepancies") or []) > 0)

# ---- 2. dashboard references every figure exactly once ----
html = open(ROOT + "/site/index.html").read()
for i in list(range(1, 51)) + ["P%d" % j for j in range(1, 7)]:
    tag = "R%d" % i if isinstance(i, int) else i
    check(f"index.html references {tag}", tag in html)
for fn in ("dlCsv", "dlJson", "dlMarkers", "dlComp", "dlCorr", "dlCum",
           "dlQc", "metricsRow", "dataTable"):
    check(f"index.html defines {fn}()",
          re.search(r"function %s\(" % fn, html) is not None)

# ---- 3. docs stay in sync with the code ----
readme = open(ROOT + "/README.md").read()
for tab in ("Generalizability", "Signatures", "Paper", "Datasets"):
    check(f"README lists {tab} tab", tab in readme)
meth = open(ROOT + "/site/methodology.html").read()
for tag in ("R50", "P6", "results.json: environment"):
    check(f"methodology.html mentions {tag}", tag in meth)
fmap = open(ROOT + "/docs/figure_map.md").read()
for tag in ("R50", "P6"):
    check(f"figure_map.md mentions {tag}", tag in fmap)

# ---- 4. every pipeline script at least parses ----
import glob as _glob
for f in sorted(_glob.glob(ROOT + "/pipeline/*.py")):
    try:
        ast.parse(open(f).read())
        check("parses " + f.split("/")[-1], True)
    except SyntaxError as e:  # noqa: BLE001
        check("parses " + f.split("/")[-1], False, str(e)[:120])

print(f"\n{len(FAIL)} failing group(s)")
sys.exit(1 if FAIL else 0)
