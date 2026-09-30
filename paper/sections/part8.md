# 8. Verification architecture

## 8.1 What the proof depends on

Figure 1 shows the flow of data and trust. Theorems 1 and 3 depend on the following, and on nothing else:

1. the mathematics of Sections 2–6 (Lemmas 4–20, Propositions 6, 10, 12, 18 and Theorem 7);
2. the exact computations that support Lemma 9, Proposition 10 and Proposition 12 (the SymPy script `exact/lemma5_exact.py` and the three scripts in `fourier/`, which use integer and rational arithmetic);
3. the checker: the kernel `chk/kernel.c` and Python modules `chk/ck.py`, `chk/run_check.py`, `chk/verify_split.py`, `chk/receipt_binding.py` and `chk/full_replay.py`, which build the family data from the defining formulas, implement the exact ordering test, drive branch and bound (including parameter subdivision), bind replay inputs, run the clique search and walk the cover tree; current line counts and hashes are in Appendix C;
4. `mpmath` interval arithmetic, used for the table verification, the constants and the coefficient rectangles of the Fourier family;
5. the usual platform: IEEE-754 binary64 hardware arithmetic with correctly rounded $+,-,\times,\sqrt{\ }$; a C compiler that honours `-ffp-contract=off` and does not reorder floating-point operations at `-O2`; the Python interpreter and NumPy for integer array manipulation and control flow (no NumPy floating-point result enters an enclosure; the family data are exact integers and the kernel receives only integers and verified binary64 end-points).

The generator, its numerical solvers, its Arb check of the table, its selector, its pair intervals and its clique enumeration are *not* in this list. The checker reads the generator's results files as plain JSON data (box indices and hull end-points); it imports no generator module. The only artefact produced by the generator that the checker uses beyond the certificates themselves is the unit table, which the checker verifies independently with a different library (`mpmath` rather than Arb).

![Data flow and trust. Theorems 1 and 3 depend only on the components in the lower (trusted) box; the generator may be wrong without affecting them, because every certificate it produces is re-verified by the checker.](fig_dependency.png){ width=100% }

## 8.2 Independence boundaries

The referee asked for a precise statement of what each verification step checks and what it shares with the main computation. Table 1 gives it.

| Question | Verified by | Shares with the generator |
|---|---|---|
| Is the unit table correct to $2^{-53}$ per component? | checker, `mpmath.iv`, 40 digits | the table file (the object being checked) |
| Does the selector keep the canonical owner of every point in $\mathcal F$? | checker's selector (Lemma 11, Section 4.3) | nothing (different formulas and arithmetic) |
| Is every kept box certified (cover completeness)? | checker's cover tree (Lemma 20) | the list of certified boxes (checked data) |
| Do the hulls contain all unbiased vectors? | checker's coverage B&B (Proposition 18) | the hull end-points (checked data) |
| Are the hulls self-separated, and is (Cl) true? | checker's pair intervals and set-based clique search | the hull end-points (checked data) |
| Are the symmetry identities and fundamental domains correct? | SymPy and exact rational scripts | nothing |

Table 1. What is verified independently of the generator.

Here “independent checker” means implementation separation from the generator. It does not mean independent human review or external reproduction of every certificate. The supplied 30 September review reports fresh graph checks for all certificates but fresh phase-space replay for seven selected parents only. Its remaining phase-space evidence was the submitted replay history.

**Replay binding (v3).** Each success now carries full SHA-256 binding to the family arrays, parameter indices and level, hull endpoints, unit table, checker and subdivision sources, build script, compiled kernel and runtime versions. Legacy or mismatched caches are rejected. `--clean` requires a new output and forbids partial-run options; proof drivers disallow `MUB6_KERNEL`. Content hashes bind inputs but do not authenticate a dishonest cache writer: the clean full replay imports no historical success records. Explicit runtime checks replace assertions on the proof path, including the family-installation call, so Python optimization cannot remove them.

The v3 regression suite retains the boundary-owner example, stale-cache driver attack, changed-input binding tests and an exact generic mixed-parameter Krawczyk fixture. The latter shows why dropping the mixed term can falsely exclude a genuine root; it is not claimed to kill the surviving mutation of the MUB kernel.

## 8.3 Tests

*Checker test-suite* (`chk/tests/test_kernel.py`, eleven groups, all passing on both machines):
T1 compares interval addition, subtraction, multiplication, squaring and square roots with exact rational arithmetic ($16{,}052$ enclosures, including the referee's value $1+2^{-52}$, zeros, subnormals and the largest finite number);
T2 checks $3{,}000$ enclosures of $e^{i\varphi}$ over angle intervals up to three full turns against 40-digit values;
T3 checks $10{,}800$ rectangles for $S_j$ at random points of random phase and parameter boxes of both families;
T4 submits $373{,}188$ dyadic phase boxes containing numerically converged solutions (residual below $10^{-10}$) to the modulus, mean-value and Krawczyk tests, none of which may exclude them;
T5 samples pair values and checks the self-separation bound and the orthogonality edges of the Fourier basis;
T6 checks that the selector keeps the boxes of $5{,}465$ (point, level) pairs of the fundamental domain, including $\alpha=0$ and both edges of the sector;
T7 compares the containment test with brute-force enumeration and with exact expectations for boxes whose end-points coincide with, or lie one grid unit beyond, those of a hull, including across the wrap-around ($950$ decisions);
T8 compares the exact ordering test with brute-force enumeration over rational sub-grids;
T9 compares the family data with the displayed formulas and checks the Hadamard property on the constraint surface;
T10 stresses the three exclusion tests with coarse parameter boxes (levels $3$–$7$) and non-dyadic phase boxes of widths $1$ to $1024$ grid units placed so that a solution lies just inside a face ($1{,}385{,}847$ boxes);
T11 checks the Krawczyk assembly directly: for sampled vertices and interior points $(w,q)$, the exactly evaluated point $(w-\tilde w)-YG(w;q)$ must lie in the kernel's interval vector $\mathcal K$ ($2{,}880$ points).

*Mutation testing* (`chk/tests/mutate.py`; Jia and Harman 2011). The supplied record reports @@MUT_KILLED@@ of @@MUT_TOTAL@@ kernel mutations detected by the test suite. The survivor removes the mixed-parameter Jacobian term in Lemma 17. This is incomplete test coverage, not evidence that the term is unnecessary or impossible to test. The exact generic fixture $G(w;q)=(1+q)w-5/2$, $W=[-1,1]$, $Q=[0,2]$, $\tilde w=0$, $c_Q=1$, $Y=1/2$ has root $(w,q)=(5/6,2)$. Omitting the mixed term gives $\mathcal K=\{5/4\}$, falsely disjoint from $W$; including it gives $[3/4,7/4]$, containing the root. This permanent rational-arithmetic test is not a claim to have killed the specific MUB-kernel mutant.

*Generator regression tests* (`tests/test_gen.py`) reproduce each counterexample of the review and check the repair: the disc product and the point–interval matrix product for $r=1+2^{-52}$ (also with zero left radius); $20{,}000$ random disc products against exact rational arithmetic; the cyclic-hull example and $3{,}000$ adversarial wrapping and nested cases; the referee's synthetic "refined" record, which is no longer counted as certified; and the selector on sampled points of $\mathcal F$ and on the special points $\alpha=0$, $\alpha=1$ and the ray.

*Historical platforms.* The submitted kernel test history records three platforms (Section 9.1): x86-64 Linux with GCC 13.3, ARM64 Linux with GCC 11.4, and ARM64 macOS with Apple clang 21. Those historical logs report all groups passing; on the two platforms where the suite was run in one piece (x86-64 Linux and macOS) every reported count is identical. They do not establish a three-platform replay of the revised Python drivers. The current local environment and fresh results are recorded separately in the release receipt. On ARM64, compilers fuse multiplications and additions by default; the explicit `-ffp-contract=off` prevents this.
