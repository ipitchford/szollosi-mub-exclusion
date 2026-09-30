# Appendix C. Inventory of code and data

Line counts and SHA-256 digests (first 16 hexadecimal digits) of the files on which the results depend; the bundle file `SHA256SUMS` lists the full digests of every file, including all results, checker outputs and logs.

**Checker (trusted)**

| file | lines | SHA-256 (prefix) |
|---|---:|---|
| `chk/kernel.c` | 638 | `9214d0d935ad21a4` |
| `chk/ck.py` | 385 | `b134fbef3bae0b90` |
| `chk/run_check.py` | 162 | `382b7910b5195adc` |
| `chk/verify_split.py` | 77 | `8ba748ea3e78101e` |
| `chk/receipt_binding.py` | 58 | `45ade521a0d9173d` |
| `chk/full_replay.py` | 55 | `556b31cae41a8769` |
| `chk/build.sh` | 11 | `fcdb7df0fe374ff1` |

**Exact verification (trusted)**

| file | lines | SHA-256 (prefix) |
|---|---:|---|
| `exact/lemma5_exact.py` | 111 | `e8142d5eb38c242f` |
| `fourier/symmetries.py` | 129 | `837092d6674b4a0a` |
| `fourier/fundamental_domain.py` | 111 | `1913323ee3ebb240` |
| `fourier/fourier_exact.py` | 17 | `6c64edd26c29f8c2` |

**Tests**

| file | lines | SHA-256 (prefix) |
|---|---:|---|
| `chk/tests/test_kernel.py` | 408 | `1d55f31c282ccb07` |
| `chk/tests/test_revision.py` | 63 | `83dfeefc5abcdd6b` |
| `chk/tests/test_boundary_interval.py` | 40 | `67470e36835c335d` |
| `chk/tests/mutate.py` | 57 | `6f504b9dfbb59629` |
| `tests/test_gen.py` | 91 | `5806d6e653e7ffec` |

**Generator (untrusted)**

| file | lines | SHA-256 (prefix) |
|---|---:|---|
| `gen/rnd.py` | 63 | `b5b23bec4e1cecff` |
| `gen/family.py` | 87 | `ceca8f8526730917` |
| `gen/core.py` | 320 | `45e90739fa9de956` |
| `gen/hulls.py` | 116 | `788317275e06297f` |
| `gen/comb.py` | 73 | `bff994877171e98f` |
| `gen/selector.py` | 104 | `eed4a8ade8ea638b` |
| `gen/solve.py` | 197 | `c692cdcf00c7fad7` |
| `gen/driver.py` | 131 | `cd4fd20bdb3215b3` |
| `gen/helpers_v2.py` | 12 | `ca1ba100a1baec90` |

**Diagnostics and analysis (not part of the proof)**

| file | lines | SHA-256 (prefix) |
|---|---:|---|
| `exact/degeneracy_diagnostics.py` | 145 | `0169911f82c81913` |
| `analysis/stats.py` | 102 | `c64b1811f6ac85d9` |
| `analysis/empty_boxes.py` | 31 | `6663ceeddaa6f7bf` |
| `analysis/ready_roots.py` | 23 | `c454da242201d947` |
| `analysis/probe_split_boxes.py` | 44 | `0f4c9d23bf6ea7cc` |
| `analysis/special_boxes.py` | 23 | `082156ef59bb789f` |

**Data**

| file | lines | SHA-256 (prefix) |
|---|---:|---|
| `gen/table17.npy` | — | `075160008692331e` |

Totals: checker (trusted), 1386 lines; exact verification (trusted), 368 lines; tests, 659 lines; generator (untrusted), 1103 lines; diagnostics and analysis (not part of the proof), 368 lines.

*Version boundary.* These hashes describe revised source. Historical logs refer to earlier versions and are not fresh receipts for these bytes. Revised clean-replay summaries bind their exact sources, compiled kernel, certificate inputs and environment. Generator variations cannot replace checker verification.

