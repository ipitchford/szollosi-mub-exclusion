# 5. The checker's arithmetic

## 5.1 One-ulp outward rounding

The checker's kernel (`chk/kernel.c`) performs all floating-point work of the proof. It uses IEEE-754 binary64 arithmetic in the default rounding mode, round-to-nearest, and never changes the mode; the only operations applied to enclosure end-points are $+,-,\times$ and $\sqrt{\ }$. It is compiled without value-changing optimisations (`-O2 -ffp-contract=off -fno-fast-math`), so that no multiply–add fusion or reassociation alters the operations written in the source. For a binary64 number $x$ let $\mathrm{up}(x)$ be the next binary64 number above $x$ and $\mathrm{dn}(x)=-\mathrm{up}(-x)$, computed by integer manipulation of the bit pattern ($\mathrm{up}(0)$ is the smallest positive subnormal, $\mathrm{up}(+\infty)=+\infty$, NaN is propagated).

**Lemma 13.** Let $a,b$ be binary64 numbers, $\circ\in\{+,-,\times\}$, and $s=\mathrm{fl}(a\circ b)$ the round-to-nearest result (possibly $\pm\infty$ after overflow). Then $\mathrm{dn}(s)\le a\circ b\le\mathrm{up}(s)$. The same holds for $s=\mathrm{fl}(\sqrt a)$, $a\ge0$.

*Proof.* If $s$ is finite, the exact result lies between the two binary64 neighbours of $s$, because round-to-nearest selects $s$ only for values at most half a spacing away from $s$ on either side, whereas the neighbours are a full spacing away (at a binade boundary the spacing below is half the spacing above, and the half-spacing criterion still refers to the correct side; in the subnormal range the spacing is constant and the argument is the same). If $s=+\infty$ the exact value exceeds the largest finite number $M$, and $\mathrm{dn}(+\infty)=M$; symmetrically for $-\infty$. $\square$

Interval operations are the standard inf–sup rules built from Lemma 13: $[a,b]+[c,d]=[\mathrm{dn}(a+c),\mathrm{up}(b+d)]$, a product takes the minimum and maximum of the four computed end-point products and steps them outward, a product containing $0\cdot\infty$ becomes $(-\infty,\infty)$, and a square is clamped at $0$ from below. Complex quantities are rectangles (pairs of intervals) with the usual rectangular product. Multiplication by an integer exponent $e$ is carried out exactly when $e\in\{0,\pm1\}$ (the term is omitted or its sign changed) and by an outward-rounded scalar product otherwise; exactly-zero structural entries of a matrix are skipped. None of these choices can make an enclosure smaller than the set it encloses.

The kernel is short enough to be read line by line (Appendix C gives its size and checksum), and its test-suite (Section 8.3) compares every operation with exact rational arithmetic on edge cases (zeros, subnormals, the largest finite number, $1+2^{-52}$) and on $4{,}000$ random pairs spread over sixty orders of magnitude.

## 5.2 The unit table and constants

All angles are integers on the grid $2\pi\mathbb Z/2^{17}$ ("units", one unit $=A:=\pi/2^{16}$). A table $T_m\approx e^{2\pi im/2^{17}}$, $0\le m<2^{17}$, supplies unit complex numbers. The table file is produced by the generator, but the checker does not rely on the generator's verification of it: before any use, the checker evaluates $\cos(2\pi m/2^{17})$ and $\sin(2\pi m/2^{17})$ for every $m$ in `mpmath` interval arithmetic at forty decimal digits and requires the rigorous upper bound of each difference from the stored component to be at most $2^{-53}$. All $262{,}144$ components pass. The kernel then stores the rectangle $[\mathrm{dn}(\mathrm{Re}\,T_m-2^{-53}),\mathrm{up}(\mathrm{Re}\,T_m+2^{-53})]\times[\cdots]$, which contains $e^{2\pi im/2^{17}}$. The worker processes verify the SHA-256 digest of the table file against the digest of the file that was verified.

Three constants enter: an upper bound $\pi^+$ of $\pi$ (the binary64 successor of $\mathrm{fl}(\pi)$), so that $A^+=\pi^+/2^{16}$ (an exact scaling) bounds $A$ from above; and binary64 numbers $\sqrt3^-<\sqrt3<\sqrt3^+$. All three inequalities are checked in $60$-digit arithmetic when the checker starts. The constant $6$ is exact.

## 5.3 Enclosures of $e^{i\varphi}$

For an integer centre $c$ and an integer half-width $w\ge0$ (units), every $\varphi\in[c-w,c+w]A$ satisfies $e^{i\varphi}=e^{icA}e^{it}$ with $|t|\le wA$, and
$$
\mathrm{Re}\,e^{it}\in[\cos wA,1]\ \ (wA\le\pi),\qquad \mathrm{Im}\,e^{it}\in[-\sin wA,\sin wA]\ \ (wA\le\pi/2),
$$
with $[-1,1]$ outside these ranges. The kernel takes $\cos wA$ and $\sin wA$ from the verified rectangle of $T_w$ (the lower end of the real part and the upper end of the imaginary part) and multiplies the rectangle of $T_{c\bmod2^{17}}$ by the rectangle $[\cos wA,1]\times[-\sin wA,\sin wA]$. The same enclosure, with $w$ the sum of the half-widths of all angles involved (the phase half-width plus $\sum_r|E_{kjr}|$ times the parameter half-width), gives a rectangle for each term
$$
\tau_{kj}(\theta,p)=\overline{c_{kj}}\,e^{i(\theta_k-\langle E_{kj},p\rangle)}
$$
of $S_j=\sum_k\tau_{kj}$ over a phase box and a parameter box, and hence a rectangle for $S_j$. Four variants are used: all angles at their centres; only the free parameters varying; the phases and the dependent parameters varying; and everything varying (Section 6.3). Derivatives are formed from the same rectangles: $\partial_{\theta_k}\tau_{kj}=i\tau_{kj}$ and $\partial_{p_r}\tau_{kj}=-iE_{kjr}\tau_{kj}$, and similarly for second derivatives and for the constraint terms $\gamma_te^{i\langle e_t,p\rangle}$.

**Lemma 14.** For every phase box $X$, parameter box $P$ and $j$, the rectangle computed for $S_j$ (in any of the four variants) contains $S_j(\theta,p)$ for all $(\theta,p)$ in the corresponding box; likewise for the derivative rectangles.

*Proof.* Lemma 13, the enclosure of $e^{i\varphi}$ above, the table verification, and the inclusion-monotonicity of interval arithmetic. $\square$

