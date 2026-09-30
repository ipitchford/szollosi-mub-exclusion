# 7. The generator

The generator (`gen/`, about $1{,}000$ lines of Python with NumPy) produces the certificates. None of its components is needed for the correctness of Theorems 1 and 3: a wrong hull family would be rejected by the checker, and a missing box would leave a hole that the checker's cover tree reports. The generator is nevertheless written to be sound in its own right, both because an unsound generator wastes time on certificates the checker rejects and because the referee's findings concerned it.

*Arithmetic.* The generator uses complex discs, NumPy products and a padding rule described in Appendix A. We withdraw the claimed general proof of that rule. The earlier counterexamples and 20,000 random disc products remain diagnostic regression tests, not universal assurance. The generator's Arb table check is likewise separate from the checker's fresh interval verification. The theorem does not depend on soundness of the generator arithmetic.

*Selector.* The generator's parameter selector implements the same three conditions as the checker's (Section 3.4) with its own formulas and disc arithmetic, and no library trigonometric function; it keeps a subset of the boxes kept by the earlier, non-rigorous selector ($14{,}658$ against $15{,}648$ level-$9$ boxes), the difference being the additional half-plane $\mathrm{Re}\,\alpha\ge0$ and the sharper remainder $6h^2$ in place of $12h^2+10^{-9}$.

*Per-box procedure.* For a root (a kept level-$7$ box) the generator runs its own modulus-only root pass (phase levels $4$ to $7$), locates the unbiased vectors at a parameter point of the root numerically (damped Newton from $20{,}000$ random starts), and then treats each kept level-$9$ child $P$:

1. the root survivors are re-tested for $P$ with a modulus test, a mean-value test and a parametric Krawczyk test (second-order Taylor form in the free parameters);
2. the numerical solutions are re-converged at a point of $P\cap\Sigma$; each receives a *declared radius* derived from a Krawczyk interval at level $L_{\max}$, clamped to $[0.02,0.15]$;
3. surviving boxes are subdivided and re-tested down to level $L_{\max}$ ($9$, and $10$ in a second attempt); from level $L_d$ ($7$, respectively $8$) on, boxes within the declared radius of a numerical solution are attached to it and not subdivided;
4. uncovered survivors are used as Newton starting points ("rescue", at most three rounds), and any still uncovered are grouped into connected clusters;
5. each solution with attached boxes, and each cluster, yields one hull: the smallest cyclic box containing all its boxes;
6. the self-separation and clique conditions are checked with the generator's own pair intervals and a bit-set clique enumeration.

Numerical solutions influence only the grouping of surviving boxes, never whether a box survives. Two defects found by the referee are repaired. The cyclic hull (step 5) now reduces every interval modulo $2^{16}$, splits wrapping intervals, merges the pieces and takes the complement of the largest gap including the wrap-around gap; the referee's example ($[1000,1100]$ and $[65000,68000]$) and $3{,}000$ adversarial cases with wrapping and nested intervals are contained correctly (comment 4.4). Hulls are formed from the attached dyadic boxes only, not from declared balls, and the Krawczyk contraction of the earlier version was dropped: on $82$ sample boxes it changed no hull size and no outcome. Records now carry an explicit state — *certified*, *refined* (the box was subdivided; its hulls are not stored) or *failed* — and on resumption only certified records are skipped (comment 4.5, with a regression test inserting the referee's synthetic record).

If a box fails at $L_{\max}=10$ it is subdivided and its kept children are treated in the same way, down to parameter level $12$. The certified records store the hulls on the grid $2\pi\mathbb Z/2^{16}$; these, and nothing else from the generator, are read by the checker.
