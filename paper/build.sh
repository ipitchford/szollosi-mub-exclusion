#!/bin/bash
# build.sh -- paper.md -> paper.tex -> paper.pdf (xelatex, Latin Modern OpenType fonts)
set -euo pipefail
cd "$(dirname "$0")"
IN=${1:-paper.md}; BASE=${IN%.md}
pandoc "$IN" -s -o "$BASE.tex" --pdf-engine=xelatex \
  -V linkcolor=blue -V urlcolor=blue -V colorlinks=true
tectonic --keep-logs "$BASE.tex"
