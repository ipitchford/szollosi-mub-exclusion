---
title: "No Szöllősi matrix lies in a quadruplet of mutually unbiased bases in dimension six: a certificate replayed by an independent checker, with a rigorous re-proof for the Fourier family"
author: "Anonymous"
date: "@@DATE@@"
geometry: "a4paper, margin=2.4cm"
fontsize: 11pt
header-includes:
  - \usepackage{amsmath,amssymb,amsthm}
  - \usepackage{booktabs,longtable}
  - \usepackage{xurl}
  - \setlength{\emergencystretch}{3em}
---

# Abstract

Matolcsi, Matszangosz, Varga and Weiner (MMVW, 2026) reduced Zauner's conjecture — that $\mathbb C^6$ carries at most three mutually unbiased bases (MUBs) — to an algebraic identity for MUB triplets (their Conjecture 2) together with their Conjecture 3: no complex Hadamard matrix in Szöllősi's two-circulant family $X(\alpha)$ is part of a quadruplet of MUBs. They remarked that Conjecture 3 should be provable with the discretisation technique of Jaming, Matolcsi, Móra, Szöllősi and Weiner (2009), "although the actual implementation could be tedious". We prove Conjecture 3 for the whole closed parameter region, boundary and degenerate loci included.

The proof is a computer-assisted certificate with a deliberately split architecture. An untrusted *generator* covers a fundamental domain of the family's symmetry group by @@SZ_BOXES@@ parameter boxes and, for each box, produces a finite family of phase-space boxes ("hulls") intended to contain every vector unbiased to both the standard basis and the Hadamard basis, for every parameter in the box. A separate *checker*, which shares no code with the generator, then replays every certificate: its trusted kernel (@@KERNEL_LINES@@ lines of C) moves every floating-point result one unit in the last place outward, uses no library transcendental function, and decides by exhaustive branch and bound (with a modulus test, a mean-value test and a parametric Krawczyk test) that the hulls cover all unbiased vectors; it then recomputes the orthogonality and unbiasedness graphs of the hulls and shows, by its own clique search, that no two further mutually unbiased bases fit into them. The checker also re-derives the parameter cover with its own selector. All @@SZ_BOXES@@ certificates were verified, @@SZ_SPLIT_COUNT@@ of them after subdividing their parameter boxes.

The same interface was applied to the two-parameter Fourier family $F(a,b)$. Exact symmetry and tiling checks justify its fundamental triangle, and @@F_BOXES@@ certificates give a new implementation of the exclusion already proved by Jaming et al. (2009), not a new Fourier exclusion theorem. The general dimension-six maximum remains conditional on MMVW's Conjecture 2. “Independent checker” refers to separation from the generator, not external peer review or proof-assistant formalisation.

# 1. Introduction

## 1.1 The problem

Two orthonormal bases $\{e_k\}$ and $\{f_l\}$ of $\mathbb C^d$ are *mutually unbiased* if $|\langle e_k,f_l\rangle|^2=1/d$ for all $k,l$. If the bases are the columns of unitary matrices $X_1,X_2$, the condition says that $\sqrt d\,X_1^*X_2$ is a *complex Hadamard matrix*: a matrix with unimodular entries whose columns are orthogonal. In prime-power dimensions there are $d+1$ mutually unbiased bases; in dimension six the largest known sets have three members, and Zauner's conjecture (Zauner 2011) asserts that three is the maximum. The question is the best-known open case of the MUB problem and has resisted analytic, numerical and algebraic attacks for more than two decades (McNulty and Weigert 2026 give a recent review).

By a unitary change of basis any set of MUBs can be brought into a form in which the first basis is the standard basis ($X_1=I$) and the second is $H/\sqrt6$ for a complex Hadamard matrix $H$. Every further basis consists of unit vectors unbiased to both. Up to a phase, the vectors unbiased to $I$ are $v(\theta)/\sqrt6$ with $v(\theta)=(1,e^{i\theta_2},\dots,e^{i\theta_6})^{\mathsf T}$, and unbiasedness to $H/\sqrt6$ becomes a system of trigonometric equations in $\theta\in\mathbb T^5$ whose solution set we denote by $\mathcal U(H)$. A quadruplet $\{I,H/\sqrt6,C,D\}$ exists exactly when $\mathcal U(H)$ contains two orthonormal bases that are unbiased to each other.

MMVW (2026) studied MUB triplets through their transition matrices and showed that one algebraic identity (their Conjecture 2) forces every Hadamard matrix occurring in a triplet to lie in one of two families: the two-parameter Fourier family $F(a,b)$, which Jaming et al. (2009) had shown cannot be extended to a quadruplet, and Szöllősi's two-parameter family $X(\alpha)$ (Szöllősi 2010). For the latter they stated:

> **Conjecture 3 (MMVW).** Any matrix $S\in X(\alpha)$ in the Szöllősi family, defined in (4.2), cannot be part of an MUB-quadruplet. In principle, this conjecture can be proven with the discretization technique in [Jaming et al. 2009], although the actual implementation could be tedious.

An earlier claim of this statement, derived from a lemma of Chen, Liang, Hu and Chen (2021), is not secure: McNulty and Weigert (2025) identified a defect in the underlying argument, and the reply of Liang and Chen (2025) does not restore it. The display (4.2) of MMVW is Szöllősi's two-circulant form
$$
S=\begin{bmatrix} A&B\\ B^*&-A^*\end{bmatrix},\qquad A=\mathrm{circ}(a,b,c),\quad B=\mathrm{circ}(d,e,f),\qquad |a|=\dots=|f|=1 ,
$$
and $X(\alpha)$ denotes the complex Hadamard matrices of this form, which up to equivalence are parametrised by a complex number $\alpha$ in a closed region $\mathbb D$ bounded by arcs of a deltoid and its negative (Section 3).

## 1.2 Results

> **Theorem 1.** Let $S$ be a complex Hadamard matrix of order six of the form (4.2). There are no unitary matrices $C,D$ such that $\{I,S/\sqrt6,C,D\}$ is a set of four mutually unbiased bases of $\mathbb C^6$. Equivalently, no member of $X(\alpha)$, for any $\alpha\in\mathbb D$ including the boundary, is part of an MUB quadruplet.

> **Corollary 2.** If Conjecture 2 of MMVW holds, then the maximal number of mutually unbiased bases in $\mathbb C^6$ is three.

The corollary is Proposition 4.2 of MMVW with Conjecture 3 discharged by Theorem 1; the Fourier half of their argument rests on Jaming et al. (2009), which we re-prove:

> **Theorem 3 (Jaming, Matolcsi, Móra, Szöllősi and Weiner 2009, Theorem 1.4).** For no $(a,b)\in[0,1)^2$ can the pair $\{I,F(a,b)\}$ be extended to a quadruplet of mutually unbiased bases.

Theorem 3 is not new. We supply a separated checker and inspectable certificates for the published Fourier exclusion. No implementation-specific criticism of the historical software is needed for this contribution. Fourier also tests transfer between the two implemented families: they share a checker and certificate format, with family-specific exponent data, coefficient rectangles and selectors. Extensions beyond these families remain future work.

## 1.3 What is new, and how it relates to Jaming et al.

The overall strategy — cover the parameter space by cells, enclose the unbiased vectors for each cell, and exclude quadruplets combinatorially — is that of Jaming et al. (2009), and MMVW explicitly anticipated its use for Conjecture 3. The contributions of this paper are the following.

1. *The theorem.* Theorem 1 is new; it removes one of the two remaining hypotheses in MMVW's route to Zauner's conjecture.
2. *A family-level certificate interface* (Section 2): an *admissible family* is specified by integer exponent vectors, unimodular coefficients and an optional trigonometric constraint; a *certificate* for a parameter box is a list of hulls in the phase torus; and four proof obligations (selector soundness, coverage, self-separation of hulls, and a clique condition) imply, for every fundamental-domain parameter in its ownership cell, that the Hadamard matrix is not part of a quadruplet. Closed boxes supply the interval enclosures; half-open cells supply ownership. The interface is what made the Fourier family a change of input data rather than a new program.
3. *Treatment of degenerate parameters without isolating roots.* Along the ray $\arg\alpha=\pi/6$ and near $\alpha=0$ unbiased vectors coalesce (Section 10). Jaming et al. separated vectors with a fixed-resolution grid; methods that isolate simple roots fail near coalescence at any fixed resolution. Here the exclusion of phase boxes is driven by a *parametric* Krawczyk operator acting uniformly over a parameter box, and hulls are allowed to contain clusters of solutions; the combinatorial step only needs the hulls to be small enough for the self-separation condition.
4. *A verification design with an explicit trusted computing base* (Sections 5–8). The generator is untrusted. The checker re-derives the owner cover, re-verifies phase-space coverage, and recomputes the graph exclusion using outward-rounded intervals. The test suite detects 16 of 17 supplied kernel mutations; that is substantial but incomplete mutation coverage. Section 8 states the surviving mutation and the additional exact generic mixed-term fixture without claiming exhaustive testing.

## 1.4 History of this work

The supplied first review reproduced counterexamples to the earlier generator's radius arithmetic, cyclic hulls and resumption logic, and identified an uncertified selector margin. The version-2 package introduced a separate checker and new certificates. The supplied second review then identified the closed-box ownership mismatch and an input-unbound replay cache. This revision repairs those contracts and withdraws the unused generator-rounding argument. Neither review demonstrated a false exclusion theorem. Review provenance is not independently authenticated, and no external human peer review is asserted. Appendix B records the second-review dispositions; historical records remain distinct from the new clean replay.

## 1.5 Organisation

Section 2 sets up admissible families, the quadruplet property and the certificate, and states the proof obligations. Sections 3 and 4 treat the Szöllősi and Fourier families: parametrisation, exact symmetries, fundamental domains and selectors. Section 5 describes the checker's arithmetic and Section 6 its exclusion tests and coverage verification. Section 7 describes the generator. Section 8 describes the verification architecture, the trusted computing base and the tests. Section 9 reports the computations and Section 10 some numerical diagnostics of degenerate parameters, which are not part of the proof. Section 11 discusses scope and limitations.
