# Appendix B. Response to the second review

This revision responds to the 30 September 2026 review of version 2. Historical logs remain evidence of the submitted version, not fresh verification of the revised driver.

| Obligation | Disposition |
|---|---|
| Closed enclosures versus half-open ownership | Theorem 7, Lemma 11 and Lemma 20 now use canonical half-open owners. Direct interval tests use closed enclosures; subdivision covers the retained fundamental-domain portion. The explicit boundary counterexample is retained as a regression. |
| Cache accepts changed hulls | Receipts bind the actual hulls, family, table, checker sources, compiled kernel and runtime. Legacy and mismatched receipts are rejected. Clean replay starts with no receipts; proof drivers reject kernel overrides. The stale-cache attack is an executable negative control. |
| Python optimization | A further defect found during revision placed the family-installation call inside an assertion. This is now an explicit runtime call; optimized-Python regression tests exercise it. |
| Unsupported larger-family transfer | Section 11.4 treats larger Karlsson and Szöllősi families as future work, with finite torus-character capacities and additional algebraic-branch obligations. |
| Invalid generator-rounding induction | Appendix A withdraws the unused lemma and explains the exponent error. Generator padding is not part of the proof; checker containment remains necessary. |
| Mixed-parameter Krawczyk term | An exact generic example demonstrates false exclusion if the term is omitted. This is not a claim to kill the surviving MUB-kernel mutation. |
| Subdivision depth | The account includes all eight historical depth-three parents and the maximum permitted level 12. |
| Historical numerical-code comparison | Unsupported implementation criticism has been removed. The Fourier precedent is credited as an existing theorem, not a new result. |
| Assurance and counts | Numerical solution counts, hull counts, graph checks and phase-space replay are distinguished. The supplied reviewer checked all graphs but only selected phase-space parents afresh. Implementation separation is not external peer review. |
| Reusable interface | A standalone certificate-interface description states ownership, arithmetic, coverage and graph obligations. Larger-family encoding and formalization remain optional future work, not achievements of this package. |

The release replay receipt, rather than historical per-machine summaries, determines the verification status of this revision.
