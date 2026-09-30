# 3. The Szöllősi family

## 3.1 Normal form and parameters

Throughout, $\mathrm{circ}(a_0,a_1,a_2)$ is the $3\times3$ circulant with rows $(a_0,a_1,a_2)$, $(a_2,a_0,a_1)$, $(a_1,a_2,a_0)$. Circulants commute, are normal, and $\mathrm{circ}(a_0,a_1,a_2)^{\mathsf T}=\mathrm{circ}(a_0,a_2,a_1)$.

**Lemma 8 (normalisation and the Hadamard condition).** (i) Multiplying rows $1$–$3$ of $S$ in (4.2) by $\lambda$ and rows $4$–$6$ by $\bar\lambda$, or columns $1$–$3$ by $\mu$ and columns $4$–$6$ by $\bar\mu$ ($|\lambda|=|\mu|=1$), gives again a matrix of the form (4.2), with $(A,B)$ replaced by $(\lambda\mu A,\lambda\bar\mu B)$; consequently every matrix of the form (4.2) is equivalent to one with $a=d=1$. (ii) A matrix of the form (4.2) with $a=d=1$ is a complex Hadamard matrix if and only if $\bar b+b\bar c+c+\bar e+e\bar f+f=0$.

*Proof.* (i) The scaled matrix is $\begin{bmatrix}\lambda\mu A&\lambda\bar\mu B\\ \bar\lambda\mu B^*&-\bar\lambda\bar\mu A^*\end{bmatrix}$, and $(\lambda\bar\mu B)^*=\bar\lambda\mu B^*$, $(\lambda\mu A)^*=\bar\lambda\bar\mu A^*$. Choosing $\lambda^2=\bar a\bar d$ and $\mu=\bar a/\lambda$ gives $\lambda\mu=\bar a$, $\lambda\bar\mu=\bar d$. (ii) Since $A,B$ are circulant, $SS^*=\mathrm{diag}(AA^*+BB^*,\,B^*B+A^*A)$ and $A^*A=AA^*$, $B^*B=BB^*$. $AA^*$ is the circulant whose first row is the autocorrelation of $(1,b,c)$, namely $(3,\ \bar b+b\bar c+c,\ \overline{\bar b+b\bar c+c})$, and likewise for $B$. Hence $SS^*=6I$ if and only if the stated sum vanishes. $\square$

For $p=(\theta_x,\theta_y,\theta_u,\theta_v)\in\mathbb T^4$ put
$$
x=e^{i\theta_x},\ y=e^{i\theta_y},\ z=\overline{xy},\qquad u=e^{i\theta_u},\ v=e^{i\theta_v},\ w=\overline{uv},
$$
$$
A(p)=\mathrm{circ}(1,\bar x,z),\quad B(p)=\mathrm{circ}(1,\bar u,w),\qquad H(p)=\begin{bmatrix}A(p)&B(p)\\ B(p)^*&-A(p)^*\end{bmatrix}.
$$
This is (4.2) with $a=d=1$, $b=\bar x$, $c=z$, $e=\bar u$, $f=w$, and every unimodular $(b,c,e,f)$ arises ($x=\bar b$, $y=b\bar c$, $u=\bar e$, $v=e\bar f$). In these variables $\bar b+b\bar c+c=x+y+z$, so by Lemma 8 $H(p)$ is Hadamard if and only if
$$
C(p):=x+y+z+u+v+w=0. \tag{C}
$$
The family is admissible in the sense of Section 2.1: every entry of $H(p)$ is $\pm e^{i\langle E,p\rangle}$ with $E\in\{-1,0,1\}^4$ having at most two non-zero entries, and $C$ has six terms. An exact symbolic computation (`exact/lemma5_exact.py`, SymPy, treating $x,y,u,v$ as Laurent variables with $\bar x=1/x$) confirms that every diagonal entry of $H(p)H(p)^*$ is $6$ and every off-diagonal entry is $0$, $C(p)$ or $\overline{C(p)}$, identically in $p$.

*Relation to Szöllősi's parameter.* Put $\alpha=x+y+z$. Since $xyz=1$ and $|x|=|y|=|z|=1$, $xy+yz+zx=\bar\alpha$, so $x,y,z$ are the roots of $f_\alpha(t)=t^3-\alpha t^2+\bar\alpha t-1$, and on $\Sigma$ the numbers $u,v,w$ are the roots of $f_{-\alpha}$. Szöllősi's region is $\mathbb D=\{\alpha: D[\alpha]\le0,\ D[-\alpha]\le0\}$ with $D[\alpha]=|\alpha|^4+18|\alpha|^2-8\,\mathrm{Re}(\alpha^3)-27$, the condition for $f_\alpha$ to have three unimodular roots (Szöllősi 2010, Section 3). Conversely, for $\alpha\in\mathbb D$ every labelling of the roots of $f_\alpha$ as $(x,y,z)$ and of $f_{-\alpha}$ as $(u,v,w)$ gives a point of $\Sigma$. Thus $\Sigma$ is precisely the set of normalised matrices of the family for all $\alpha$ in the *closed* region $\mathbb D$ and all labellings; points on the boundary of $\mathbb D$, where roots coincide, are points of $\Sigma$ like any others and need no separate treatment.

## 3.2 Symmetries

Let $\Pi=\mathrm{circ}(0,1,0)$, let $P_{23}$ exchange indices $2$ and $3$, $D_\omega=\mathrm{diag}(1,\omega,\omega^2)$ with $\omega=e^{2\pi i/3}$, and write $M_1\oplus M_2$ for a block-diagonal matrix. Define the following maps of $\mathbb T^4$, written through the six unimodular numbers $(x,y,z,u,v,w)$, which determine $p$:
$$
\begin{aligned}
\sigma_1&:(x,y,z,u,v,w)\mapsto(y,z,x,u,v,w), & \sigma_2&:(x,y,z,u,v,w)\mapsto(x,y,z,v,w,u),\\
\rho&:(x,y,z,u,v,w)\mapsto(\bar z,\bar y,\bar x,\bar w,\bar v,\bar u), & \tau&:(x,y,z,u,v,w)\mapsto(\bar z,\bar y,\bar x,\bar u,\bar v,\bar w),\\
\kappa&:(x,y,z,u,v,w)\mapsto(\bar x,\bar y,\bar z,\bar u,\bar v,\bar w), & \omega_\cdot&:(x,y,z,u,v,w)\mapsto\omega\,(x,y,z,u,v,w),\\
\varsigma&:(x,y,z,u,v,w)\mapsto(u,v,w,x,y,z). &&
\end{aligned}
$$

**Lemma 9 (matrix identities).** For every $p\in\mathbb T^4$:

1. $H(\tau p)=H(p)^{\mathsf T}$;
2. $H(\kappa p)=\overline{H(p)}$;
3. $H(\rho p)=(P_{23}\oplus P_{23})\,H(p)\,(P_{23}\oplus P_{23})$;
4. $H(\omega_\cdot p)=(D_\omega\oplus D_\omega)\,H(p)\,(D_\omega\oplus D_\omega)^{-1}$;
5. $H(\varsigma p)=(I\oplus(-I))\,H(p)\,\begin{bmatrix}0&I\\I&0\end{bmatrix}$;
6. $H(\sigma_1p)=(xI\oplus\Pi)\,H(p)\,(\Pi^{\mathsf T}\oplus\bar xI)$ and $H(\sigma_2p)=(uI\oplus\Pi)\,H(p)\,(\bar uI\oplus\Pi^{\mathsf T})$.

Moreover $C(\sigma_1p)=C(\sigma_2p)=C(\varsigma p)=C(p)$, $C(\rho p)=C(\tau p)=C(\kappa p)=\overline{C(p)}$ and $C(\omega_\cdot p)=\omega C(p)$, so every map preserves $\Sigma$; and $\alpha$ is mapped to $\alpha$ by $\sigma_1,\sigma_2$, to $\bar\alpha$ by $\rho,\tau,\kappa$, to $\omega\alpha$ by $\omega_\cdot$, and to $u+v+w$ ($=-\alpha$ on $\Sigma$) by $\varsigma$.

*Proof.* Write $A=\mathrm{circ}(1,\bar x,z)$, $B=\mathrm{circ}(1,\bar u,w)$. (1) $H^{\mathsf T}=\begin{bmatrix}A^{\mathsf T}&\bar B\\ B^{\mathsf T}&-\bar A\end{bmatrix}$ with $A^{\mathsf T}=\mathrm{circ}(1,z,\bar x)$ and $\bar B=\mathrm{circ}(1,u,\bar w)$; this pair is $(\mathrm{circ}(1,\bar x',z'),\mathrm{circ}(1,\bar u',w'))$ with $x'=\bar z$, $z'=\bar x$, $y'=\bar y$, $u'=\bar u$, $w'=\bar w$, $v'=\bar v$, and the lower blocks are the required adjoints. (2) is immediate. (3) Conjugation by $P_{23}$ transposes a circulant; the rest is as in (1). (4) $D_\omega\,\mathrm{circ}(a_0,a_1,a_2)\,D_\omega^{-1}=\mathrm{circ}(a_0,\bar\omega a_1,\bar\omega^2a_2)$. (5) $H\begin{bmatrix}0&I\\I&0\end{bmatrix}=\begin{bmatrix}B&A\\-A^*&B^*\end{bmatrix}$, and negating the last three rows gives $H(\varsigma p)$. (6) $xA=\mathrm{circ}(x,1,xz)=\mathrm{circ}(x,1,\bar y)$ is a cyclic shift of $A(\sigma_1p)=\mathrm{circ}(1,\bar y,x)$; the shifts commute with circulants and cancel on the lower-right block. All identities, including the exponent signs in (6), which the earlier version left as "a suitable choice", and the statements about $C$ and $\alpha$, were verified by exact symbolic computation (`exact/lemma5_exact.py`, output `lemma5_exact.json`). $\square$

By Lemma 4, $\mathcal E(H(gp))\iff\mathcal E(H(p))$ for every $g$ in the group $G$ generated by the seven maps.

## 3.3 A fundamental domain

Use angle coordinates $\theta_x,\dots,\theta_w\in[0,2\pi)$ for the six numbers, with $\theta_z\equiv-\theta_x-\theta_y$ and $\theta_w\equiv-\theta_u-\theta_v$. Consider

* (S) $\alpha=x+y+z$ lies in the closed sector $\{0\le\arg\alpha\le\pi/6\}\cup\{0\}$;
* (OA) $\theta_x\le\theta_y\le\theta_z$; (OB) $\theta_u\le\theta_v\le\theta_w$.

**Proposition 10.** For every $p\in\Sigma$ there is $g\in G$ such that $gp$ satisfies (S), (OA) and (OB). Hence if $\mathcal E(H(p))$ fails for every $p$ in $\mathcal F:=\{p\in\mathbb T^4:\ (\mathrm C),(\mathrm S),(\mathrm{OA}),(\mathrm{OB})\}$, it fails for every $p\in\Sigma$.

*Proof.* $\sigma_1$ and $\sigma_2$ permute the labels of the two triples cyclically. The composition $\kappa\tau$ maps $(x,y,z,u,v,w)$ to $(z,y,x,u,v,w)$ and $\kappa\rho$ maps it to $(z,y,x,w,v,u)$ (both verified symbolically), so $(\kappa\rho)(\kappa\tau)$ reverses the second triple alone. Cyclic permutations and a reversal generate $S_3$; hence $G$ contains $S_3\times S_3$ acting on the labels, which leaves $\{x,y,z\}$, $\{u,v,w\}$ and $\alpha$ unchanged, and sorting the angles of each triple achieves (OA) and (OB). The maps $\omega_\cdot$, $\varsigma$ and $\kappa$ send $\alpha$ to $\omega\alpha$, $-\alpha$ and $\bar\alpha$; they generate the dihedral group of order twelve (rotations by multiples of $\pi/3$ and a reflection), whose closed fundamental sector is $0\le\arg\alpha\le\pi/6$ ($\alpha=0$ is fixed). After moving $\alpha$ into the sector, a relabelling, which does not change $\alpha$, restores (OA) and (OB). All maps preserve (C). $\square$

## 3.4 The selector

Ownership cells and their closed enclosures are as in Section 2.1. The enclosure indexed by $i$ at level $\ell$ has centre $c=(2i+1)\pi/2^\ell$ and half-width $h=\pi/2^\ell$. Centres lie on $2\pi\mathbb Z/2^{17}$ for $\ell\le16$. The selector rejects only when a test excludes the ownership cell from $\mathcal F$; the ordering test need not exclude its entire closed enclosure.

*Ordering (OA), (OB): exact integer test.* In units of $2\pi/2^\ell$ write $t_1\in[i_1,i_1+1)$, $t_2\in[i_2,i_2+1)$ for $\theta_x,\theta_y$. Every point of $\mathcal F$ has one half-open owner at each level. Direct interval tests use its closed enclosure; when subdivision is needed, only the retained fundamental-domain portion is certified by following the unique ownership path. Then $\theta_z$ corresponds to $t_3=(-t_1-t_2)\bmod N$, $N=2^\ell$. (OA) requires $t_1\le t_2\le t_3$, hence $i_1\le i_2$ and $\sup t_3\ge i_2$. The quantity $-(t_1+t_2)$ ranges over $(\beta-2,\beta]$ with $\beta=-(i_1+i_2)$, which contains the integers $\beta-1$ and $\beta$; if one of them is a multiple of $N$ (that is, $\beta\bmod N\in\{0,1\}$) the residue $t_3$ approaches $N$ and $\sup t_3=N$, otherwise $\sup t_3=\beta\bmod N$. The box is discarded if $i_1>i_2$ or $\sup t_3<i_2$; (OB) is treated identically. The generator uses an algebraically different formula for the same condition; the checker's test-suite compares its formula with brute-force enumeration on fine rational sub-grids.

*Constraint (C): support-function test with rigorous rounding.* Let $D(c)=(i(x-z),\,i(y-z),\,i(u-w),\,i(v-w))$ be the gradient of $C$ at the centre (complex numbers viewed in $\mathbb R^2$). For $|\delta_r|\le h$ Taylor's theorem gives $C(c+\delta)=C(c)+\sum_rD_r(c)\delta_r+R$ with $|R|\le6h^2$: the Hessian of each term $e^{i\langle e_t,p\rangle}$ in direction $\delta$ has modulus at most $\langle e_t,\delta\rangle^2\le\|e_t\|_1^2h^2$, and $\sum_t\|e_t\|_1^2=1+1+4+1+1+4=12$, so the second-order remainder is at most $\tfrac12\cdot12h^2$. If $C(c+\delta)=0$ then for every vector $n\in\mathbb R^2$
$$
\langle n,-C(c)\rangle\le h\sum_r|\langle n,D_r(c)\rangle|+6h^2|n| .
$$
The box is discarded if some $n$ violates this inequality. The candidate vectors $n$ (edge normals and vertex directions of the zonotope $\{\sum_rD_r\delta_r\}$ and $-C(c)$ itself, $25$ in all) are computed in ordinary floating point, which is harmless because the inequality holds for *every* $n$; what must be rigorous is its evaluation, and in the checker every term is evaluated in the interval arithmetic of Section 5: $C(c)$ and $D_r(c)$ are rectangles built from verified table entries, $|n|$ is bounded above with outward steps, $h\ge\pi/2^\ell$ is an upward-rounded product with an upper bound of $\pi$, and the comparison is between a rigorous lower bound of the left side and a rigorous upper bound of the right side.

*Sector (S).* Over the box, $|x-x(c)|,|y-y(c)|\le h$ and $|z-z(c)|\le2h$ (chord at most arc), so $\alpha$ lies within $4h$ of $\alpha(c)$, which lies in a rigorously computed rectangle. The closed sector is contained in each of the half-planes $\{\mathrm{Im}\,\alpha\ge0\}$, $\{\mathrm{Re}\,\alpha\ge0\}$ and $\{\mathrm{Re}\,\alpha-\sqrt3\,\mathrm{Im}\,\alpha\ge0\}$. The box is discarded if the disc of radius $4h$ about the rectangle lies strictly outside one of them: $\max\mathrm{Im}+4h<0$, or $\max\mathrm{Re}+4h<0$, or $\min(\sqrt3\,\mathrm{Im}-\mathrm{Re})>2\cdot4h$ (the normal $(-1,\sqrt3)$ has length $2$), with $\sqrt3$ replaced by the neighbouring binary64 numbers below and above it (verified with 60-digit arithmetic) according to the sign of $\mathrm{Im}$.

**Lemma 11 (owner retention).** The selector retains the canonical owner of every $p\in\mathcal F$ at each relevant level. A discarded half-open ownership cell does not meet $\mathcal F$; its closed enclosure may do so.

*Proof.* The ordering inequalities above are necessary for the half-open owner of a point satisfying (OA)/(OB). The constraint and sector tests exclude only when their condition fails throughout the closed enclosure, hence throughout the ownership cell. A point in $\mathcal F$ fails none of these necessary conditions in its owner. $\square$

At level 9, take $t=\pi/32$, $u=v=-e^{it}$, $w=e^{-2it}$, $\alpha=2e^{it}-e^{-2it}$ and ordered roots $x,y,z$ of $f_\alpha$. The resulting fundamental-domain point has retained owner $(34,106,264,264)$. It also lies on the upper $v$ boundary of the discarded closed enclosure $(34,106,264,263)$. This permanent regression illustrates why closed-enclosure retention is not asserted.

These tests replace the selector of the earlier version, which evaluated centres with the library exponential and relied on a margin of $10^{-9}$ (referee comment 4.2). Starting from level $3$, the checker's selector keeps @@SZ_KEPT_LEVELS@@ boxes at levels $3,4,\dots,9$. The level-$9$ boxes are the ones certified; the @@SZ_ROOTS@@ kept level-$7$ boxes ("roots") organise the computation (Section 7).
