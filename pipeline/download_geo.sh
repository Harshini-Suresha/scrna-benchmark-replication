#!/usr/bin/env bash
# GSE133549 is a SuperSeries without files of its own: data live in 15 SubSeries (one per protocol).
# Usage: bash pipeline/download_geo.sh list   |   bash pipeline/download_geo.sh get GSE133540 [GSE...]
set -euo pipefail
ALL="GSE133535 GSE133536 GSE133537 GSE133538 GSE133539 GSE133540 GSE133541 GSE133542 GSE133543 GSE133544 GSE133545 GSE133546 GSE133547 GSE133548 GSE141469"
url() { echo "https://ftp.ncbi.nlm.nih.gov/geo/series/${1:0:6}nnn/$1/suppl/"; }
ls_files() { curl -fsS "$(url $1)" | grep -o 'href="[^"?/][^"]*"' | cut -d'"' -f2 || true; }
case "${1:-list}" in
  list) for a in $ALL; do echo "== $a"; ls_files $a; done ;;
  get)  shift; for a in "${@:-$ALL}"; do mkdir -p data/raw/$a
        for f in $(ls_files $a); do case "$f" in http*|*.html|*.php) continue;; esac; curl -fSL -o data/raw/$a/$f "$(url $a)$f"; done; done ;;
esac
echo "If a SubSeries lists no supplementary files, use its Series Matrix or the SRA Run Selector (PRJNA551748)."
