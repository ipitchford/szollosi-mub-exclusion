# Szöllősi-family MUB exclusion

**Unrefereed candidate, revision 3.0.** Fresh local replay passed all 14,658
Szöllősi boxes and all 1,935 Fourier boxes, including the required boundary
subdivisions. See [STATUS.md](STATUS.md) and [replay receipt](REPLAY_RECEIPT.md).
This is not external peer review or end-to-end formal verification.

The candidate excludes a quadruplet of mutually unbiased bases in complex
dimension six containing the standard basis and a matrix from the specified
closed two-circulant Szöllősi family, including boundary and degenerate points.
It does **not** settle the unrestricted dimension-six MUB problem. The latter
consequence remains conditional on MMVW's separate Conjecture 2. The package's
Fourier exclusion re-proves the established theorem of Jaming et al. (2009).

## Start here

- [Manuscript PDF](paper/paper.pdf), [Markdown](paper/paper.md), [TeX](paper/paper.tex).
- [AI index](AI_INDEX.md): claims, evidence, replay and reuse boundaries.
- [Assurance](ASSURANCE.md), [certificate interface](CERTIFICATE_INTERFACE.md),
  [source audit](CITATION_AUDIT.md), [provenance](PROVENANCE.md).
- Appendix B of the manuscript answers the second review. Earlier response
  files and README-v2-HISTORICAL.md are superseded where they conflict with revision 3.

## Reproduce

Use Python 3.10 or later, NumPy, mpmath, SymPy and a C99 compiler. The published
environment record specifies the local versions actually used. The generator
additionally uses python-flint; it is not needed for certificate replay.

```sh
sh -e chk/build.sh
python3 chk/tests/test_kernel.py
python3 chk/tests/test_revision.py
python3 -O chk/tests/test_revision.py
python3 -O chk/tests/test_boundary_interval.py
python3 exact/lemma5_exact.py
python3 fourier/fundamental_domain.py
python3 fourier/fourier_exact.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 chk/full_replay.py /path/to/new-output 4
```

The output directory must not already exist. Full replay imports no historical
success receipts, checks both families, subdivides unresolved parents, and
requires both retained parameter covers to be complete. Budget hours rather
than minutes for the full phase-space calculation. A regenerated kernel has a
new build identity; do not mix its receipts with those of another build.

`sh -e` makes compilation fail closed. Do not set `MUB6_KERNEL` for proof runs.
Do not run `final_eval.sh` to validate revision 3: it is retained only as a
disabled historical driver for the submitted computation.

## Evidence layout

| Path | Role |
|---|---|
| `chk/` | Trusted interval kernel, selectors, cover/graph checks and receipt binding. |
| `exact/`, `fourier/` | Exact symmetry and Fourier-domain checks. |
| `gen/`, `results/` | Untrusted generator and proposed certificates; no root-count assurance. |
| `chk_results/`, `logs/` | Submitted historical replay records, not revised-driver acceptance. |
| `replay/` | Fresh input-bound revision-3 replay, complete summaries, subdivisions and original binary. |
| `analysis/` | Historical statistics and non-proof diagnostics. |
| `paper/sections/`, `paper/values.json` | Revised manuscript sources and explicitly historical numerical statistics. |
| `review/` | Supplied earlier review material, with separate licence/provenance boundary. |

The repaired argument follows canonical half-open parameter owners while
performing interval tests on closed enclosures. Receipts bind actual family
arrays, parameter boxes, hulls, table, sources, binary and runtime. Hashes detect
mismatches, not fraudulent receipt authors; use clean replay for reproduction.
The unused generator-rounding lemma is withdrawn. Tests detect 16 of 17 supplied
kernel mutations, not every possible defect. Implementation separation from the
generator is not external peer review or formal verification.

## Rebuild the manuscript

```sh
python3 paper/inventory.py
python3 paper/assemble.py
bash paper/build.sh
```

Pandoc and Tectonic are required. The assembly does not import old replay passes.
Preserve the historical statistics unless explicitly replacing their labelled
history with a new measured run.

## Reuse

Scholarly creator: Anonymous. Original prose, diagrams and data: CC0-1.0.
Original code: MIT. Third-party and reviewer material retains its own terms;
see [LICENSES.md](LICENSES.md). No licence conveys a correctness warranty.
