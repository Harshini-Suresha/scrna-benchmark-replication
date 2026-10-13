set -e
python pipeline/make_synthetic.py && python pipeline/run_analysis.py config.synthetic.yaml
python - <<'PY'
import json; d=json.load(open("site/data/demo.json")); assert len(d["protocols"])==3
for p in d["protocols"]: assert p["ari"]>0.5, p
print("OK")
PY
