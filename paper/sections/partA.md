# Appendix A. Withdrawn generator rounding claim

The generator is outside the trusted computing base. The former Lemma R and its proof are withdrawn: multiplication combines uncertainty from both operands, invalidating the claimed intermediate induction. With $u=2^{-53}$, two exact inputs $(1-u)^{-8}$ represented by computed values $1$ have product $(1-u)^{-16}$, exceeding the asserted exponent $9$ even with its small additive allowance. This does not refute the final padded bound, but that unused bound is no longer claimed. The padding rule below remains a generator heuristic; every proposal must be checked independently.

The implemented rule is $\mathrm{up}(x)=\mathrm{fl}(\mathrm{fl}(x(1+2^{-40}))+2^{-800})$.

Generator regression tests are diagnostics, not proof of a general rounding lemma. The theorem relies on the separate outward-rounded checker, not these heuristics.

Every radius returned by a generator primitive is passed through $\mathrm{up}$ (and replaced by $+\infty$ if it exceeds $2^{64}$ or is not a number); differences and lower bounds are formed by a single correctly rounded operation on exact operands followed by a one-ulp step, as in Lemma 13; and all exclusion decisions are strict comparisons of a lower bound with an upper bound, which are false whenever a NaN or infinity occurs.
