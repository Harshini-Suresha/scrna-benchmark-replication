"""Resolve marker-gene symbols via mygene.info (graceful offline fallback).

Reads site/data/results.json, collects union marker ENSEMBL IDs (versions
stripped for lookup), batch-resolves symbols, and stores
protocol.markers.symbols = {as_stored_id: symbol}. Dashboard prefers symbols,
falls back to IDs. Rewrites results.js too. Safe to re-run (idempotent).
Usage: python pipeline/map_symbols.py [results.json]
"""
import json
import sys
import urllib.request

RES = sys.argv[1] if len(sys.argv) > 1 else "site/data/results.json"
CHUNK = 200


def resolve(ids):
    out = {}
    ids = list(ids)
    for i in range(0, len(ids), CHUNK):
        batch = ids[i : i + CHUNK]
        req = urllib.request.Request(
            "https://mygene.info/v3/gene",
            data=json.dumps({"ids": ",".join(batch), "fields": "symbol"}).encode(),
            headers={"Content-Type": "application/json"},
        )
        try:
            for hit in json.load(urllib.request.urlopen(req, timeout=60)):
                if "symbol" in hit and "query" in hit:
                    out[hit["query"]] = hit["symbol"]
        except Exception as e:
            print(f"  batch {i // CHUNK + 1} failed ({e}); continuing with IDs")
    return out


def main():
    d = json.load(open(RES))
    want = set()
    for p in d["protocols"]:
        mk = p.get("markers") or {}
        for g in mk.get("genes", []):
            want.add(g.split(".")[0])
    print(f"{len(want)} unique marker IDs to resolve")
    if not want:
        print("no markers; nothing to do")
        return
    syms = resolve(sorted(want))
    print(f"resolved {len(syms)}/{len(want)}")
    for p in d["protocols"]:
        mk = p.get("markers") or {}
        if mk.get("genes"):
            mk["symbols"] = {
                g: syms[g.split(".")[0]]
                for g in mk["genes"]
                if g.split(".")[0] in syms
            }
    json.dump(d, open(RES, "w"))
    bn = RES; var = "window.DEMO" if "demo" in bn else ("window.RESULTS_MOUSE" if "mouse" in bn else "window.RESULTS")
    open(RES.rsplit(".", 1)[0] + ".js", "w").write(
        var + " = " + json.dumps(d, separators=(",", ":")) + ";"
    )
    print("updated", RES)


main()
