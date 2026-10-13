"""Maps downloaded SubSeries folders (data/raw/<GSE>/) to <outdir>/<protocol>.<ext>.
Env: PROCESSED_DIR (default data/processed), SPECIES=human|mouse (default human; filters to "*_human_*" / "*_mouse_*" tables, falling back to the largest table).
Extracts tar archives, picks the LARGEST table file per SubSeries and prints the alternatives; verify with inspect_data.py."""
import glob, os, shutil, sys, tarfile
MAP = {"GSE133535": "10X_2x5K_250K", "GSE133536": "10X_8x10K_25K", "GSE141469": "10X_V3", "GSE133537": "C1HT_medium",
       "GSE133538": "C1HT_small", "GSE133539": "CELseq2", "GSE133540": "Dropseq", "GSE133541": "ICELL8", "GSE133542": "MARSseq",
       "GSE133543": "QUARTZseq2", "GSE133544": "gmSCRBseq", "GSE133545": "SMARTseq2", "GSE133546": "SingleNuclei",
       "GSE133547": "ddSEQ", "GSE133548": "inDrop"}
EXT = (".txt.gz", ".tsv.gz", ".csv.gz", ".txt", ".tsv", ".csv")
OUTDIR = os.environ.get("PROCESSED_DIR", "data/processed")
SPECIES = os.environ.get("SPECIES", "human")
os.makedirs(OUTDIR, exist_ok=True)
for acc in sys.argv[1:] or MAP:
    d = f"data/raw/{acc}"
    for t in glob.glob(d + "/*.tar*"):
        tarfile.open(t).extractall(d)
    c = sorted([f for f in glob.glob(d + "/**/*", recursive=True) if f.endswith(EXT)], key=os.path.getsize, reverse=True)
    if not c: print(f"{acc}: no table files in {d}"); continue
    sp = [f for f in c if f"_{SPECIES}_" in os.path.basename(f)]
    if sp: c = sp
    ext = next(e for e in EXT if c[0].endswith(e))
    print(f"{acc} -> {MAP[acc]}{ext}   chosen: {os.path.basename(c[0])}   others: {[os.path.basename(x) for x in c[1:4]]}")
    shutil.copy(c[0], f"{OUTDIR}/{MAP[acc]}{ext}")
