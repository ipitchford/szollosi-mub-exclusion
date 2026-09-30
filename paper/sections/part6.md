# 6. Exclusion tests and the verification of a certificate

Throughout, $X$ is a phase box $\prod_{k=2}^6[l_k,h_k]\cdot2\pi/2^{16}$ (centre $l_k+h_k$ and half-width $h_k-l_k$ in units of $A$), $P$ a closed parameter box with centre $c_P$ and half-width $h_P$, and "$X\times P$ contains no solution" means: there is no $(\theta,p)\in X\times P$ with $p\in\Sigma$ and $\theta\in\mathcal U(H(p))$. Only necessary conditions for being a solution are used to exclude, so each test is valid whether or not $H(p)$ is Hadamard at the points of $P$ off $\Sigma$.

## 6.1 Modulus test

For each $j$ the kernel forms the rectangle $M_j$ of $\sum_k\tau_{kj}$ with all angles at their centres, and the radius $R_j=A^+\sum_kw_{kj}$ with $w_{kj}$ the total half-width of the angle of $\tau_{kj}$ (upward rounded). Since $|e^{i\varphi}-e^{i\psi}|\le|\varphi-\psi|$, $S_j(\theta,p)$ lies within $R_j$ of $M_j$.

**Lemma 15.** If for some $j$ the upper bound of $|M_j|+R_j$ has square below $6$, or the lower bound of $\mathrm{dist}(0,M_j)-R_j$ is positive with square above $6$, then $X\times P$ contains no solution. $\square$

This test is cheap and discards most of the phase torus; it is also the only test used in the *root pass*, which subdivides the whole phase torus uniformly from level $4$ ($2^{20}$ boxes) to level $7$ for a level-$7$ parameter box and keeps the survivors (between @@ROOT_MIN@@ and @@ROOT_MAX@@ boxes at level $7$ for the Szöllősi roots). A box excluded for a parameter box is excluded for every sub-box, so the root pass is shared by all certified boxes below the same root.

## 6.2 Mean-value test

Let $F_j=|S_j|^2-6$. The box $X\times P$ is convex, so by the mean-value theorem
$$
F_j(\theta,p)\in F_j(\tilde\theta,c_P)+\sum_{k=2}^6\partial_{\theta_k}F_j(X\times P)\,[-\rho_k,\rho_k]+\sum_{r}\partial_{p_r}F_j(X\times P)\,[-h_P,h_P],
$$
where $\tilde\theta$ is the centre of $X$, $\rho_k$ the half-widths of $X$ in radians (upward rounded), and $\partial F_j=2\,\mathrm{Re}(\overline{S_j}\,\partial S_j)$ is enclosed over the whole box using Lemma 14.

**Lemma 16.** If for some $j$ the interval above does not contain $0$, then $X\times P$ contains no solution. $\square$

## 6.3 Parametric Krawczyk test

For a family with a constraint the parameters are split, at the centre of $P$, into $n_d=2$ *dependent* angles $d$ and the remaining *free* angles $q$; any split is valid, and the checker takes the pair for which the $2\times2$ block of the differential of $C$ at $c_P$ has the largest determinant (a floating-point choice that affects only efficiency). For a family without a constraint $n_d=0$ and all parameters are free. The unknowns are $w=(\theta_2,\dots,\theta_6,d)\in W=X\times D$, $n=5+n_d$, and the equations are
$$
G(w;q)=\big(F_2,\dots,F_6,\ (\mathrm{Re}\,C,\ \mathrm{Im}\,C)\big)\in\mathbb R^{n}
$$
(the last two only if there is a constraint). Every solution $(\theta,p)$ with $p\in P\cap\Sigma$ is a zero $w^*\in W$ of $G(\cdot;q)$ for the corresponding $q\in Q$. Let $\tilde w$ be the centre of $W$, $c_Q$ that of $Q$, $\rho$ the (upward-rounded) half-widths of $W$, and $Y$ any real $n\times n$ matrix; the kernel takes a floating-point inverse of the Jacobian at $(\tilde w;c_Q)$ (and makes no exclusion if the inversion fails).

**Lemma 17.** Let $F_q(\tilde w;Q)$ be an interval matrix containing $\partial_qG(\tilde w;q)$ for all $q\in Q$, $J(W;c_Q)$ one containing $\partial_wG(w;c_Q)$ for all $w\in W$, and $J_q(W;Q)$ one containing $\partial_q\partial_wG(w;q)$ for all $(w,q)\in W\times Q$. Put
$$
V=Y\,G(\tilde w;c_Q)+\sum_l\big(Y F_q(\tilde w;Q)\big)_{\cdot l}\,[-h_P,h_P],
$$
$$
\mathcal M=I-Y J(W;c_Q)-\sum_l\big(YJ_{q_l}(W;Q)\big)[-h_P,h_P],
$$
$$
\mathcal K=-V+\mathcal M\,[-\rho,\rho].
$$
If for some $r$ the interval $\mathcal K_r$ is disjoint from $[-\rho_r,\rho_r]$, then $X\times P$ contains no solution.

*Proof.* Fix $q\in Q$ and suppose $G(w^*;q)=0$ with $w^*\in W$. Applying the mean-value theorem to each component of $g(w)=w-YG(w;q)$ on the segment from $\tilde w$ to $w^*$ gives $w^*-\tilde w=g(w^*)-\tilde w=-YG(\tilde w;q)+(I-\Lambda)(w^*-\tilde w)$, where the $r$-th row of $\Lambda$ is $\sum_sY_{rs}\,\partial_wG_s(\xi_r;q)$ for some point $\xi_r$ of the segment, which lies in $W$. By the mean-value theorem in $q$, applied entrywise, $G(\tilde w;q)\in G(\tilde w;c_Q)+F_q(\tilde w;Q)(q-c_Q)$ and $\partial_wG(\xi;q)\in J(W;c_Q)+J_q(W;Q)(q-c_Q)$ for $\xi\in W$, with $|q_l-c_{Q,l}|\le h_P$. Because $Y$ is a point matrix, $Y\big(F_q(q-c_Q)\big)\subseteq\sum_l(YF_q)_{\cdot l}(q_l-c_{Q,l})$, and likewise $\Lambda\in YJ(W;c_Q)+\sum_l(YJ_{q_l})(q_l-c_{Q,l})$ entrywise; hence $w^*-\tilde w\in\mathcal K$. Also $w^*-\tilde w\in[-\rho,\rho]$. So $\mathcal K_r\cap[-\rho_r,\rho_r]\ne\emptyset$ for every $r$. $\square$

The kernel evaluates all enclosures with Lemma 14 and all products and sums in outward-rounded interval arithmetic, applying $Y$ before multiplying by $[-h_P,h_P]$ so that the correlation between the equations captured by $Y$ is retained. The earlier version used a second-order Taylor form in $q$ for the value term; the first-order mean-value form above is simpler to verify and was sufficient.

## 6.4 Containment, coverage and the pair graphs

*Containment.* A box $[l,h]$ (per coordinate, grid units) lies in a hull $[\ell,\eta]$ if, for every coordinate with $\eta-\ell<2^{16}$, the offset $o=(l-\ell)\bmod2^{16}\in[0,2^{16})$ satisfies $o+(h-l)\le\eta-\ell$. This integer test can only err by declaring a box not contained.

*Coverage verification.* For a certified parameter box $P$ with hull family $\mathcal H_P$ and root $P_7\supseteq P$, the checker proceeds as follows.

1. Take the survivors of the root pass for $P_7$ at level $7$ and apply the modulus test for $P$.
2. At the current level, discard every box contained in a hull of $\mathcal H_P$; apply the mean-value test and then the Krawczyk test to the others and discard the boxes they exclude.
3. If no box remains, report success. Otherwise, if the level is $12$ or more than $3\cdot10^5$ boxes remain, report failure; else replace every remaining box by its $32$ children, discard those excluded by the modulus test for $P$, and return to step 2.

**Proposition 18.** If the coverage verification reports success for $P$, then $\mathcal U(H(p))\subseteq\bigcup\mathcal H_P$ for every $p\in P\cap\Sigma$ (for every $p\in P$ if there is no constraint).

*Proof.* Let $p\in P\cap\Sigma$ and $\theta\in\mathcal U(H(p))$. The root pass starts from the closed dyadic boxes of level $4$, which cover the torus, and removes only boxes excluded for $P_7\supseteq P$ (Lemma 15), so $\theta$ lies in a surviving box, and at each subdivision in one of its children. A box is removed in step 2 only if it lies in a hull, or if it contains no solution for any parameter of $P\cap\Sigma$ (Lemmas 15–17). Hence, as long as $\theta$ is not in a hull, it lies in some box of the current list; when the list is empty, $\theta$ lies in a hull. $\square$

*Pair intervals and graphs.* For hulls $\mathsf h,\mathsf h'$ with centres $c,c'$ and half-widths $n,n'$ (units), $\langle v(\theta),v(\phi)\rangle=1+\sum_{k=2}^6e^{i(\phi_k-\theta_k)}$ with $\phi_k-\theta_k\in[c'_k-c_k\mp(n_k+n'_k)]A$; the kernel encloses this sum in a rectangle $\mathcal T$ (Section 5.3) and $|\langle v,w\rangle|^2$ in the interval $\Pi=|\mathcal T|^2$. The edge $\mathsf h\sim_O\mathsf h'$ is set if $\inf\Pi\le0$ and $\mathsf h\sim_B\mathsf h'$ if $6\in\Pi$; the self-separation (Sep) is certified if, with half-widths $2n_k$ about the centre $0$, $\inf\mathrm{Re}\,\mathcal T>0$ and $\inf|\mathcal T|^2>6$.

**Lemma 19.** The graphs so constructed contain the graphs $O$, $B$ of Section 2.4 for the true pair intervals, and a hull certified by the self-separation test satisfies (Sep). $\square$

Extra edges only make (Cl) harder to satisfy, so Lemma 19 is the direction needed. The checker enumerates all six-cliques of $O$ by a set-based recursive search and, for each, searches for a six-clique of $O$ inside the common $B$-neighbourhood of its vertices; (Cl) holds if none is found. Neither the enumeration nor the graphs share code with the generator, which uses bit-set enumeration on graphs built by different formulas (Section 7).

## 6.5 The cover tree

**Lemma 20.** Let the checker's selector (Sections 3.4 and 4.3) be applied from level $3$: the kept boxes of level $3$ form the first generation, and a kept box that is not certified is replaced by its kept children. If this process terminates, by level $12$, with every remaining box certified and verified (Proposition 18 together with the checks of Lemma 19 and the clique search), then (O1)–(O4) of Theorem 7 hold.

*Proof.* Every $p\in\mathcal F$ has one half-open owner at each level, retained by Lemma 11 or Section 4.3. Its next owner is a child of its present owner, including grid boundaries and the torus seam. The path reaches a certified cell. Proposition 18 on its enclosure implies (O2) on the ownership cell; the remaining checks give (O3)–(O4). $\square$

*Verification by parameter subdivision.* Obligation (O2) is needed only for $p\in\mathcal F\cap Q$, where $Q$ is the parent ownership cell and $P=\overline Q$. Each such point follows a retained descendant owner within $Q$. When coverage cannot be certified on all of $P$, the checker tests the same hull family on the closed enclosures of all retained descendants at levels $\ell+1$, $\ell+2$, or $\ell+3$, up to level 12. Proposition 18 on these enclosures proves coverage on $\mathcal F\cap Q$. Hull-only (Sep) and (Cl) then establish (O2)–(O4). Boundary points owned by another parent are not claimed to follow this parent's descendants. Eight submitted certificates required $\ell+3$.
