# Szöllősi family and MUB quadruplets: certificate, independent checker and evidence

**Revision in preparation, 30 September 2026.** The revised clean full replay
is still running. Historical receipts listed below are retained provenance,
not passes for the revised driver. The Fourier clean replay and normal and
optimized revision tests have passed locally; publication awaits the complete
Szöllősi replay. This is an unrefereed candidate, not externally verified.
See `CERTIFICATE_INTERFACE.md` for the current clean replay command and trust
boundary. Appendix A withdraws the old generator rounding lemma. Canonical
half-open ownership replaces the incorrect closed-box selector assertion.

This bundle accompanies the paper *No Szöllősi matrix lies in a quadruplet of mutually unbiased bases in dimension six: a certificate replayed by an independent checker, with a rigorous re-proof for the Fourier family* (`paper/paper.pdf`; Markdown source `paper/paper.md`, LaTeX `paper/paper.tex`). The response to the review of the earlier version is `paper/response_to_review.pdf` (source `.md`).

## What is proved, and where the evidence is

Theorem 1 of the paper states that no complex Hadamard matrix of the Szöllősi two-circulant form, for any parameter in the closed region including its boundary and degenerate points, is part of a quadruplet of mutually unbiased bases (MUBs) in dimension six. Theorem 3 is a rigorous re-proof of the corresponding result of Jaming, Matolcsi, Móra, Szöllősi and Weiner (2009) for the two-parameter Fourier family. Corollary 2 records that, together with Conjecture 2 of Matolcsi, Matszangosz, Varga and Weiner, this settles the question of the maximal number of MUBs in dimension six.

The computational part of each proof is the statement that the independent checker verified every certificate and that its cover evaluation is complete. The files that record this are:

| Claim | File |
|---|---|
| Szöllősi cover complete; 14,658 certificates verified (75 of them by parameter subdivision) | `chk_results/chk_szollosi_summary.json` (field `cover_complete`), produced by the run logged in `logs/chk_szollosi_final.log` |
| Fourier cover complete; 1,935 certificates verified | `chk_results/chk_fourier_summary.json`, `logs/chk_fourier.log` |
| Exact symmetry identities of the Szöllősi family (Lemma 9) | `exact/lemma5_exact.py`, output `exact/lemma5_exact.json` (the file name predates the final numbering) |
| Fourier symmetry group of order 144 and the fundamental triangle (Proposition 12) | `fourier/symmetries.py`, `fourier/fundamental_domain.py`, output `fourier/fundamental_domain.json`; exact Hadamard check `fourier/fourier_exact.py` |
| Kernel test-suite (T1–T11) on three platforms | `chk/tests/test_kernel.py`, logs in `chk/tests/platform_logs/` |
| Mutation testing of the kernel (16 of 17 mutations detected) | `chk/tests/mutate.py`, `chk/tests/mutation_results.json`, `chk/tests/mutation.log` |
| Regression tests for every counterexample of the review | `tests/test_gen.py` |

Every number quoted in the paper is computed from these outputs by `final_eval.sh` (statistics in `analysis/`, values in `paper/values.json`) and substituted into the text by `paper/fill.py`; nothing is typed in by hand.

## Layout

| Directory | Contents | Role |
|---|---|---|
| `paper/` | paper (Markdown, LaTeX, PDF), Figure 1 and its script, response to the review; `sections/` holds the section templates with `@@…@@` placeholders, `values.json` the values substituted | — |
| `chk/` | the **independent checker**: `kernel.c` (interval kernel with one-ulp outward rounding), `ck.py`, `run_check.py`, `verify_split.py` (verification by parameter subdivision), `build.sh`; `tests/` holds the kernel test-suite, the mutation tests and their logs | trusted |
| `exact/` | exact SymPy verification of the Szöllősi symmetry identities; numerical diagnostics of degenerate parameters (`degeneracy_diagnostics.py`, not part of the proof) | trusted (symmetries only) |
| `fourier/` | exact symmetry search and fundamental-domain proof for the Fourier family; exact Hadamard check | trusted |
| `gen/` | the generator: family data, arithmetic of Lemma R, selector, per-box certificate, driver; the unit table `table17.npy` | untrusted |
| `tests/` | generator regression tests; an early trial run of the generator | — |
| `results/` | raw generator output, one JSON record per parameter box: `gen_vm.jsonl` and `gen_vm_tail.jsonl` (machine A), `gen_mac.jsonl` (machine B), `gen_fourier.jsonl` (machine A), `fourier_trial.jsonl` (trial) | data checked by the checker |
| `chk_results/` | checker output, one line per root (`chk_szollosi_mac.jsonl`, `chk_szollosi_vmrec_mac.jsonl`, `chk_fourier.jsonl`), the parameter-subdivision verifications (`chk_szollosi_split*.jsonl`, with their box lists `split_boxes_*.json`), the merged file `chk_szollosi_all.jsonl` and the summaries `*_summary.json`; `superseded_split_runs/` holds three early subdivision runs that are not used | evidence |
| `analysis/` | statistics scripts and their JSON output, including the comparison of the rigorous selector with the earlier one (`selector_diff.json`) and the sampling probe of the subdivided boxes | not part of the proof |
| `logs/` | every run log (generator, checker, subdivision runs, final evaluation); `superseded_split_runs/` holds the logs of the unused early subdivision runs | evidence |
| `scripts_machineB/`, `pipeline_vm.sh`, `pipeline_vm2.sh` | the job scripts actually run on the two machines (with their original absolute paths) | provenance |
| `review/` | the referee report and the referee's audit bundle of the earlier version | — |
| `SHA256SUMS` | SHA-256 checksum of every file in the bundle | provenance |

## Reproducing the verification

Requirements: Python 3.10 or later, NumPy, `mpmath`, SymPy and a C99 compiler. The generator additionally uses `python-flint` (Arb) for its own check of the unit table; the checker does not. Set `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` when running several processes.

```
cd chk && ./build.sh && cd ..
python3 chk/tests/test_kernel.py           # kernel test-suite T1-T11
python3 chk/tests/mutate.py                # mutation tests (runs the suite once per mutation)
python3 tests/test_gen.py                  # generator regression tests

# exact computations
python3 exact/lemma5_exact.py
python3 fourier/symmetries.py; python3 fourier/fundamental_domain.py; python3 fourier/fourier_exact.py

# replay every certificate (about 15 CPU-hours for the Szöllősi family, 0.4 for Fourier)
G="results/gen_vm.jsonl results/gen_vm_tail.jsonl results/gen_mac.jsonl"
python3 chk/run_check.py szollosi gen/table17.npy chk_results/replay.jsonl 4 $G
python3 chk/run_check.py fourier  gen/table17.npy chk_results/replay_f.jsonl 4 results/gen_fourier.jsonl

# boxes whose certificate the checker cannot verify on the whole parameter box
# are verified on their kept sub-boxes (Section 6.5 of the paper); the 75 boxes are listed in split_boxes_*.json
python3 chk/verify_split.py szollosi gen/table17.npy chk_results/replay.jsonl $G --boxes chk_results/split_boxes_vmrec.json --depth 3
#   ... likewise for split_boxes_mac.json and split_boxes_mac2_0.json ... split_boxes_mac2_5.json

# cover evaluation and all statistics of the paper, from the checker outputs in the bundle
./final_eval.sh
```

`run_check.py` resumes from its output file and evaluates the cover tree at the end of every run; a later successful record for a box supersedes an earlier failed one, so appending the subdivision records to the replay file and running `run_check.py` again on it (with the same generator files) reproduces the final cover evaluation. The checker verifies the unit table with `mpmath` interval arithmetic before use, re-derives the parameter cover with its own selector, and records in its summary the SHA-256 digests of its own code and of the table.

To rebuild the paper from the templates: `cd paper && cat sections/part1.md … sections/partS.md > paper_template.md && python3 fill.py values.json paper_template.md paper.md && ./build.sh paper.md` (the section order is 1–11, A, B, C, S; `build.sh` needs pandoc and XeLaTeX).
