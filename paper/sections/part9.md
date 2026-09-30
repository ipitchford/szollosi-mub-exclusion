# 9. The computations

This section preserves the submitted generation and replay history and its measured statistics. These historical machine runs are distinct from the revised clean replay described in the release receipt. Historical cache-based summaries are not by themselves verification of the revised driver; the new release requires fresh input-bound coverage checks.

## 9.1 Set-up

Three machines were used. *Machine A* is a cloud virtual machine (two Intel Xeon cores at 2.1 GHz, 8 GB, Ubuntu 24.04, Python 3.11, NumPy 2.4.4 with OpenBLAS, `mpmath` 1.3.0, `python-flint` 0.9.0, GCC 13.3). *Machine B* is an Apple M2 Mac mini (eight cores, 16 GB, macOS, Python 3.14, NumPy 2.4.2, `mpmath` 1.3.0, `python-flint` 0.9.0, Apple clang 21), shared with other work during the run. *Machine C* is an ARM64 Linux virtual machine on the same Mac (Ubuntu 22.04, Python 3.10, NumPy 2.2.6, GCC 11.4), used for the test-suite. BLAS threading was disabled in every process. The unit table used by all runs has SHA-256 `@@TABLE_SHA@@`; it was verified by the generator with Arb and by the checker with `mpmath`, independently, on each machine.

The Szöllősi roots (kept level-$7$ boxes) were ordered by decreasing $|\alpha|$ at their centres. Machine A generated the certificates for the roots at positions $\equiv6\pmod7$ (@@SZ_VM_ROOTS@@ roots, two processes); machine B generated those at the other positions (six processes), except the last $250$ of its list, which machine A generated after finishing its share. The Fourier certificates were generated on machine A. The checker was run on machines A and B, each replaying part of the certificates, and the final cover evaluation was run on the union of all checker outputs.

## 9.2 The Szöllősi family

The checker's selector keeps @@SZ_KEPT_LEVELS@@ parameter boxes at levels $3,\dots,9$; the @@SZ_BOXES@@ kept level-$9$ boxes lie below @@SZ_ROOTS_WITH_CHILDREN@@ of the @@SZ_ROOTS@@ roots (the other roots have no kept child). @@SZ_GEN_OUTCOME@@

*Historical checker replay.* The submitted checker history reports @@SZ_CHK_VERIFIED@@. Its cover evaluation, run on the union of its outputs, reports every kept box of level $3$ covered (`chk_results/chk_szollosi_summary.json`). For the revised release, this historical summary is not the acceptance gate: Theorem 1 additionally requires the fresh input-bound full replay recorded in `REPLAY_RECEIPT.md`, together with the semantic arguments of Theorem 7, Lemma 20 and Proposition 10.

@@SZ_SPLIT_TEXT@@

Table 2 summarises the certificates. The number of numerically located solutions, which only steers the grouping of boxes into hulls and is not a certified count, is given as a complete histogram.

| quantity | value |
|---|---|
@@SZ_TABLE@@

Table 2. Statistics of the Szöllősi certificates (generator) and of their replay (checker).

The checker's coverage branch and bound needed phase levels up to @@SZ_DEEPEST@@. The total generator time was @@SZ_GEN_HOURS@@ CPU-hours and the checker time @@SZ_CHK_HOURS@@ CPU-hours.

*Special parameters.* @@SZ_SPECIAL@@

## 9.3 The Fourier family

@@F_TEXT@@

| quantity | value |
|---|---|
@@F_TABLE@@

Table 3. Statistics of the Fourier certificates and their historical replay.

## 9.4 History of the run

Three incidents are recorded because the raw results files contain their traces; none of them touches the checker, and every record written in these circumstances is either a certificate, which the checker verifies like any other, or a failed or refined record, which is not used.

First, machine B's Python installation lacked SciPy, which the generator used for the connected components that group left-over boxes into extra hulls. The resulting exceptions were treated, correctly, as failed certificates, and the affected level-$9$ boxes were subdivided down to level $12$, where @@SZ_FAILED@@ records ended in the state *failed*. The component routine was replaced by a SciPy-free union–find in NumPy (tested against SciPy on $300$ random graphs), the run was resumed, and every affected level-$9$ box was then certified at level $9$; the checker's cover tree stops at those certificates.

Second, machine B's generator was restarted twice: once with a limit on its list of roots, to hand the last $250$ roots to machine A, and once for the repair just described. The resumption logic skipped every box already certified. On macOS the worker processes of an interrupted run are not terminated with their parent, and those of one interrupted run continued for about half an hour, appending records, before they were found and stopped.

Third, machine A's container was restarted by its host during the tail run, which then resumed from its results file.

## 9.5 Relation to the earlier computation

The earlier computation certified $15{,}648$ level-$9$ boxes selected by the non-rigorous selector; the rigorous selector keeps $14{,}658$, all of them among the earlier ones. None of the earlier certificates is used. As a test of the checker's strength on the hardest parameters, it was also applied to six earlier certificates near $\alpha=0$, on the ray and near $\alpha=1$, and verified all of them.
