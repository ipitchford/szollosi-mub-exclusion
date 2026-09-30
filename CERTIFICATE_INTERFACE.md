# Certificate interface

This interface describes the two implemented families, not arbitrary algebraic
Hadamard families. See Theorem 7 and Sections 3–6 of the manuscript.

1. **Family semantics.** Supply the torus-character exponent and coefficient
   arrays defining the matrix and its constraint equations. Verify separately
   that these arrays describe the claimed family and that symmetries reduce
   its parameters to the stated fundamental domain. Current kernel capacities
   are MAXP=4, MAXT=8 and MAXN=7; algebraic branches are not silently supported.
2. **Ownership.** Dyadic parameter cells are half open with canonical periodic
   representatives. Interval arithmetic uses their closed enclosures. Every
   fundamental-domain point must retain its unique owner at each subdivision.
   Discarding an owner does not imply its closed enclosure misses the domain.
3. **Vector coverage.** Each parameter certificate supplies hull endpoints for
   candidate unbiased vectors. The generator is untrusted. Interval branch and
   bound must prove coverage of all relevant vectors; failure or exhaustion is
   not a pass. Subdivision must cover the retained fundamental-domain portion.
4. **Graph obstruction.** Verified pair enclosures establish possible
   orthogonality edges. Self-separation and the specified clique obstruction
   must hold. A graph check alone does not prove vector coverage.
5. **Global closure.** The checker cover tree must have no uncovered retained
   owner. For a split parent, all necessary retained descendants must be
   certified. Boundary and degenerate parameters are not removed by sampling.
6. **Replay identity.** Receipts bind family arrays, parameter cells, hulls,
   table, source, binary and runtime. These hashes detect mismatch, not a
   malicious receipt author. Use a new clean output directory to reproduce
   without trusting historical passes. Proof drivers reject MUB6_KERNEL.

## Reproduce

```sh
sh -e chk/build.sh
python3 chk/tests/test_revision.py
python3 -O chk/tests/test_revision.py
python3 chk/tests/test_kernel.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 chk/full_replay.py /path/to/new-output 4
```

Exit zero from the full replay requires both family covers, including any
needed subdivisions. Preserve the entire output directory and summaries.
Neither this interface nor a schema validation is a formal proof assistant
verification, independent human review, or proof of the unrestricted MUB(6)
conjecture.
