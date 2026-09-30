# Local clean replay receipt

Recorded: 2026-09-30T06:26:43.079301+00:00. Producer-side execution on one local machine.

The revised clean driver was run into a new output directory. No supplied
historical success record was imported. Required boundary subdivisions were
computed afresh and then combined with this same run for final cover checking.

```sh
sh -e chk/build.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 chk/full_replay.py NEW_OUTPUT_DIRECTORY 4
```

| Family | Required boxes | Verified | Cover |
|---|---:|---:|---|
| fourier | 1935 | 1935 | complete |
| szollosi | 14658 | 14658 | complete |

Szöllősi boundary parents replayed by subdivision: 75.
Successful subdivision levels: 38 parents at level 10, 29 parents at level 11, 8 parents at level 12.
Both final summaries have zero failed, unchecked or missing required leaves.

## Bound evidence

- `replay/fourier_summary.json` and `replay/szollosi_summary.json`: full source, table, certificate and environment identity.
- `replay/*.jsonl`: per-root/per-box input-bound results; historical failed parent attempts remain visible before successful subdivisions.
- `replay/*-base.log`, `replay/*-split-*.log`, and the final Szöllősi log: observed computation records.
- `replay/compiled/libmubk.dylib`: exact macOS binary used for this run, not a portable executable.

Kernel SHA-256: `089fceb601e22c869dbbd9611e70cbb35e02cffb98be4a2906766618e6d38450`.
Table SHA-256: `075160008692331ecec57ebd851514642c5800b39d9bcad11c668c01074c5969`.
Runtime: 3.14.7 (main, Aug  5 2026, 10:29:49) [Clang 21.0.0 (clang-2100.1.1.101)]; NumPy 2.4.2; mpmath 1.3.0; macOS-27.0-arm64-arm-64bit-Mach-O.

## Assurance ceiling

This is local replay, not external independent reproduction, human peer
review or end-to-end formal verification. Written ownership, symmetry and
encoding arguments, the interval kernel, compiler and hardware remain
separate trusted layers. Hashes bind bytes; they do not authenticate an
arbitrary fabricated pass record. `check_receipt_integrity.py` checks
recorded identity and summary completeness only; to reproduce the actual
mathematical computation run the clean driver above.

The manuscript distinguishes the new Szöllősi-family candidate from the
Fourier theorem already proved in 2009. The unrestricted dimension-six
three-basis conclusion remains conditional on MMVW Conjecture 2.
