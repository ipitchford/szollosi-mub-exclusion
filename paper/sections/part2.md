# 2. Families, certificates and proof obligations

## 2.1 Admissible families

Write $\mathbb T=\mathbb R/2\pi\mathbb Z$. An *admissible family* of order six with $m$ parameters consists of

* integer exponent vectors $E_{kj}\in\mathbb Z^m$ and unimodular constants $c_{kj}$ ($1\le k,j\le6$), defining
$$
H(p)_{kj}=c_{kj}\,e^{i\langle E_{kj},p\rangle},\qquad p\in\mathbb T^m ;
$$
* optionally, a *constraint* $C(p)=\sum_{t}\gamma_t\,e^{i\langle e_t,p\rangle}$ with integer vectors $e_t$ and constants $\gamma_t$, whose zero set $\Sigma=\{p: C(p)=0\}$ is the set of parameters at which $H(p)$ is a complex Hadamard matrix (without a constraint, $\Sigma=\mathbb T^m$ and $H(p)$ must be Hadamard for all $p$);
* a *fundamental domain* $\mathcal F\subseteq\Sigma$ with the property that the quadruplet question for every $p\in\Sigma$ reduces to the same question for some point of $\mathcal F$ (Sections 3.3 and 4.3), together with a *selector* that retains the canonical half-open dyadic owner of every point of $\mathcal F$ at every relevant level.

Use representatives in $[0,2\pi)^m$ and define
$$Q_{\ell,i}=\prod_r[2\pi i_r/2^\ell,2\pi(i_r+1)/2^\ell),\qquad P_{\ell,i}=\overline{Q_{\ell,i}}.$$
Closure is on the torus; the endpoint $2\pi$ is represented by $0$ for ownership. Each point has a unique owner at each level, and its finer owner is a child of its coarser owner. Closed enclosures are used for interval calculations. A discarded enclosure may meet $\mathcal F$ on a shared boundary.

For the Szöllősi family $m=4$, all $c_{kj}=\pm1$, and the constraint has six terms (Section 3). For the Fourier family $m=2$, the $c_{kj}$ are sixth roots of unity, and there is no constraint (Section 4). The checker consumes exactly this data: the exponent arrays, rectangles enclosing the constants, and the selector.

## 2.2 The quadruplet property

For a complex Hadamard matrix $H$ of order six let $\mathcal E(H)$ be the statement: *there exist unitary $C,D$ such that $\{I,H/\sqrt6,C,D\}$ is a set of four mutually unbiased bases.* Two Hadamard matrices are *equivalent* if $H'=M_1HM_2$ with monomial unitary matrices $M_1,M_2$ (products of a permutation matrix and a diagonal unitary matrix).

**Lemma 4.** (i) If $H'$ is equivalent to $H$ then $\mathcal E(H)\iff\mathcal E(H')$. (ii) $\mathcal E(H)\iff\mathcal E(\overline H)$. (iii) $\mathcal E(H)\iff\mathcal E(H^*)$, hence $\mathcal E(H)\iff\mathcal E(H^{\mathsf T})$.

*Proof.* (i) Let $\{I,M_1HM_2/\sqrt6,C,D\}$ be a quadruplet. Multiplying all four unitaries on the left by $M_1^{-1}$ preserves mutual unbiasedness and gives $\{M_1^{-1},HM_2/\sqrt6,M_1^{-1}C,M_1^{-1}D\}$. The columns of the monomial unitary $M_1^{-1}$ are the standard basis vectors up to phases, and the columns of $HM_2$ are those of $H$ up to order and phases; as bases these are $I$ and $H/\sqrt6$. The converse is the same argument with $M_1^{-1},M_2^{-1}$. (ii) Conjugate all four matrices; conjugation is an involution. (iii) Put $U=H/\sqrt6$ and let $\{I,U,C,D\}$ be a quadruplet. Then $\{I,U^*,U^*C,U^*D\}$ is one: $U^*$, $U^*C$, $U^*D$ have entries of modulus $1/\sqrt6$ because $U$ is unbiased to $I$, $C$ and $D$; $(U^*)^*(U^*C)=C$ and $(U^*)^*(U^*D)=D$ are unbiased to $I$; and $(U^*C)^*(U^*D)=C^*D$ is unbiased because $C$ and $D$ are. Since $(H^*)^*=H$ this is an equivalence, and $H^{\mathsf T}=\overline{H^*}$ with (ii) gives the last claim. $\square$

## 2.3 Unbiased vectors

For a complex Hadamard matrix $H$ and $\theta=(\theta_2,\dots,\theta_6)\in\mathbb T^5$ put $v(\theta)=(1,e^{i\theta_2},\dots,e^{i\theta_6})^{\mathsf T}$ and
$$
S_j(\theta)=(H^*v(\theta))_j=\sum_{k=1}^6\overline{H_{kj}}\,v_k(\theta),\qquad j=1,\dots,6 .
$$
The unit vector $v(\theta)/\sqrt6$ is unbiased to the columns of $H/\sqrt6$ if and only if $|S_j(\theta)|^2=6$ for all $j$. Let
$$
\mathcal U(H)=\{\theta\in\mathbb T^5:\ |S_j(\theta)|^2=6,\ j=1,\dots,6\}.
$$
(For a Hadamard matrix $\sum_j|S_j|^2=\|H^*v\|^2=36$, so one of the six equations is implied by the others; the tests of Section 6 use either all six or five of them, which is legitimate because a necessary condition for membership in $\mathcal U(H)$ is used to exclude.) No finiteness of $\mathcal U(H)$ is assumed anywhere.

**Lemma 5.** If $\mathcal E(H)$ holds, there are twelve distinct points $\theta^{(1)},\dots,\theta^{(12)}\in\mathcal U(H)$ such that $\langle v(\theta^{(a)}),v(\theta^{(b)})\rangle=0$ for $a\ne b$ both in $\{1,\dots,6\}$ or both in $\{7,\dots,12\}$, and $|\langle v(\theta^{(a)}),v(\theta^{(b)})\rangle|^2=6$ for $a\le6<b$.

*Proof.* Each column of $C$ and of $D$ is unbiased to $I$, so all its entries have modulus $1/\sqrt6$. Multiplying a column by a unimodular phase changes neither orthogonality within its basis nor unbiasedness to other bases; after such a phase and the factor $\sqrt6$ the column is $v(\theta)$ with $\theta\in\mathcal U(H)$. Orthogonality within $C$ and $D$ gives the zero inner products, and unbiasedness of $C$ and $D$, $|\langle c_a,d_b\rangle|^2=1/6$ for unit vectors, becomes $|\langle v,w\rangle|^2=6$ after rescaling by $\sqrt6$ each. Two of the twelve points coincide only if $|\langle v,v\rangle|^2=36$ equals $0$ or $6$, which is impossible. $\square$

## 2.4 Hulls and the combinatorial certificate

A *hull* is a box $\prod_{k=2}^{6}[l_k,h_k]\cdot\frac{2\pi}{2^{16}}$ in the phase torus with integer end-points $l_k\le h_k$ (if $h_k-l_k\ge2^{16}$ the factor is the whole circle). For hulls $\mathsf h,\mathsf h'$ let $\Pi(\mathsf h,\mathsf h')\subseteq[0,\infty)$ be an interval containing $|\langle v(\theta),v(\phi)\rangle|^2$ for all $\theta\in\mathsf h$, $\phi\in\mathsf h'$, and $\Pi(\mathsf h,\mathsf h)$ one containing these values for all $\theta,\phi\in\mathsf h$. Define two graphs on a finite family $\mathcal H$ of hulls:
$$
\mathsf h\sim_O\mathsf h'\iff \mathsf h\ne\mathsf h',\ 0\in\Pi(\mathsf h,\mathsf h'),\qquad \mathsf h\sim_B\mathsf h'\iff \mathsf h\ne\mathsf h',\ 6\in\Pi(\mathsf h,\mathsf h').
$$

**Proposition 6 (combinatorial certificate).** Let $H$ be a complex Hadamard matrix and $\mathcal H$ a finite family of hulls such that

* (Cov) $\mathcal U(H)\subseteq\bigcup_{\mathsf h\in\mathcal H}\mathsf h$;
* (Sep) $\Pi(\mathsf h,\mathsf h)\subset(6,\infty)$ for every $\mathsf h\in\mathcal H$;
* (Cl) there are no two six-cliques $K_1,K_2$ of the graph $O$ such that every pair in $K_1\times K_2$ is an edge of $B$.

Then $\mathcal E(H)$ fails.

*Proof.* Suppose $\mathcal E(H)$ holds, take the twelve points of Lemma 5 and assign to each a hull containing it (Cov). For two distinct points of the twelve, $|\langle v,w\rangle|^2\in\{0,6\}$; by (Sep) they cannot share a hull, so the twelve hulls are distinct. The six hulls of the first basis are pairwise joined in $O$, likewise those of the second basis, and every cross pair is joined in $B$ (by the definition of $\Pi$). This contradicts (Cl). $\square$

The certificate is *uniform in the parameter*: if (Cov) holds for $H(p)$ for every $p$ in a set $P'$, and (Sep) and (Cl) hold for the hull family (they do not involve $p$), then $\mathcal E(H(p))$ fails for every $p\in P'$.

## 2.5 Proof obligations

A *certificate* for a closed dyadic parameter box $P$ is a finite hull family $\mathcal H_P$. A *cover* is a set of certified boxes. The following statement is what the checker establishes.

**Theorem 7.** Let an admissible family with fundamental domain $\mathcal F$ and selector be given, and a cover with certificates $\mathcal H_P$. Suppose

* (O1) *selector soundness and cover completeness:* every kept cell of the selector's cover tree (Section 6.5) is either certified or has all its kept children covered, down to a finite depth; every $p\in\mathcal F$ has its canonical owner retained at each relevant level;
* (O2) *coverage:* for every certified ownership cell $Q$ with enclosure $P$ and every $p\in\mathcal F\cap Q$, $\mathcal U(H(p))\subseteq\bigcup\mathcal H_P$;
* (O3) *separation:* (Sep) holds for every hull of every $\mathcal H_P$;
* (O4) *clique condition:* (Cl) holds for every $\mathcal H_P$.

Then $\mathcal E(H(p))$ fails for every $p\in\mathcal F$, and hence, by the defining property of $\mathcal F$, for every $p\in\Sigma$.

*Proof.* Follow the unique owners of $p\in\mathcal F$ from the top level. By (O1) they are retained, nested, and reach a certified cell $Q$ at finite depth. Thus $p\in\mathcal F\cap Q$, and (O2)–(O4) give Proposition 6 with $\mathcal H=\mathcal H_P$. $\square$

In the Szöllősi case the reduction from $\Sigma$ to $\mathcal F$ is Proposition 10, and every normalised matrix of the form (4.2) is $H(p)$ for some $p\in\Sigma$ (Lemma 8); together these give Theorem 1. In the Fourier case the reduction is Proposition 12, which gives Theorem 3.
