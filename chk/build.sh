#!/bin/sh
# Build the checker kernel.  No -ffast-math; floating-point contraction (FMA fusion) disabled; the code never
# changes the rounding mode, and every endpoint is moved one ulp outward explicitly (see kernel.c header).
cd "$(dirname "$0")"
CC=${CC:-cc}
case "$(uname)" in
  Darwin) OUT=libmubk.dylib; SH="-dynamiclib" ;;
  *)      OUT=libmubk.so;    SH="-shared" ;;
esac
$CC -O2 -fPIC $SH -std=c99 -Wall -Wno-misleading-indentation -ffp-contract=off -fno-fast-math -frounding-math -o $OUT kernel.c -lm
echo "built $OUT with $($CC --version | head -1)"
