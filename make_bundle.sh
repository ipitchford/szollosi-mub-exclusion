#!/bin/bash
# make_bundle.sh -- assemble the evidence bundle (a single zip) from the working tree.
set -e
SRC="$(cd "$(dirname "$0")" && pwd)"
DEST=${1:-$SRC/../bundle_v2}
NAME=mub6_szollosi_certified_v2
UPLOADS=/root/.claude/uploads/057db7aa-29f4-5207-b661-48717aedb21a
OUT=$DEST/$NAME
rm -rf "$DEST"; mkdir -p "$OUT"
cd "$SRC"
# paper and response
mkdir -p $OUT/paper/sections
cp paper/paper.md paper/paper.tex paper/paper.pdf paper/fig_dependency.png paper/fig_dependency.pdf $OUT/paper/
cp paper/make_figure.py paper/build.sh paper/fill.py paper/make_values.py paper/make_values2.py paper/inventory.py \
   paper/values_base.json paper/values.json $OUT/paper/
cp paper/response_filled.md $OUT/paper/response_to_review.md
cp paper/response_filled.pdf $OUT/paper/response_to_review.pdf
cp paper/part*.md $OUT/paper/sections/
cp paper/response_to_review.md $OUT/paper/sections/response_template.md
# checker
mkdir -p $OUT/chk/tests
cp chk/kernel.c chk/ck.py chk/run_check.py chk/verify_split.py chk/build.sh $OUT/chk/
cp chk/tests/test_kernel.py chk/tests/mutate.py chk/tests/*.log chk/tests/mutation_results.json $OUT/chk/tests/
cp -r chk/tests/platform_logs $OUT/chk/tests/
# generator, exact computations, tests, analysis
mkdir -p $OUT/gen $OUT/exact $OUT/fourier $OUT/tests $OUT/analysis
cp gen/*.py gen/table17.npy $OUT/gen/
cp exact/*.py exact/*.json exact/*.log $OUT/exact/
cp fourier/*.py fourier/*.json $OUT/fourier/
cp tests/* $OUT/tests/
cp analysis/*.py analysis/*.json analysis/*.log $OUT/analysis/
# raw results, checker outputs, logs, job scripts
cp -r results chk_results logs scripts_machineB $OUT/
cp README.md final_eval.sh make_bundle.sh pipeline_vm.sh pipeline_vm2.sh $OUT/
# the review of the earlier version
mkdir -p $OUT/review
cp $UPLOADS/ff8d8b52-attachment.txt $OUT/review/referee_report.md
cp $UPLOADS/465c58aa-mub6_referee_audit_2026-09-29.zip $OUT/review/mub6_referee_audit_2026-09-29.zip
# clean, checksum, zip
find $OUT -name "__pycache__" -prune -exec rm -rf {} +
find $OUT \( -name "*.so" -o -name "*.dylib" -o -name "*.pyc" -o -name ".DS_Store" \) -delete
(cd $OUT && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
(cd $DEST && zip -qr -9 $NAME.zip $NAME)
ls -la $DEST/$NAME.zip
