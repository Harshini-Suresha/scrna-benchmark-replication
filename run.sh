#!/usr/bin/env bash
# One entry point. Usage: bash run.sh setup|download|prepare|inspect|analyse|test|serve
set -euo pipefail
case "${1:-help}" in
  setup)    python3 -m venv .venv && . .venv/bin/activate && pip install -q -r requirements.txt && echo "Activate with: source .venv/bin/activate" ;;
  download) shift; bash pipeline/download_geo.sh get "$@" ;;
  prepare)  python pipeline/prepare_data.py "${@:2}" ;;
  inspect)  python pipeline/inspect_data.py ;;
  analyse)  python pipeline/run_analysis.py ;;
  test)     bash tests/test_pipeline.sh ;;
  serve)    echo "http://localhost:8000  (add ?demo for synthetic preview)"; python -m http.server -d site 8000 ;;
  *) sed -n '2p' "$0" ;;
esac
