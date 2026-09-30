# Review: Szöllősi-family MUB exclusion

## 1. Executive assessment

**Recommendation: major revisions. Prospective REF assessment: 2★ as submitted; a credible 3★ candidate after successful repair and reverification. A 4★ assessment is not presently justified.**

This is a substantial computer-assisted mathematics manuscript with a coherent proof architecture and unusually useful supporting material. Its proposed result—excluding the entire Szöllősi family from quadruplets of mutually unbiased bases, including boundary and degenerate parameters—is clearly stated. The strongest feature is the separation between numerical solution-finding, rigorous enclosure, and a finite combinatorial obstruction: numerical roots guide the computation but need not constitute an exhaustive solution list.

My checks provide substantial positive evidence. All supplied checksums matched; every entry of the transcendental table passed an independent interval check; all **6,628,140 stored six-cliques** were re-enumerated without discrepancies; and the **15,648 final records** exactly matched the boxes retained by the supplied parameter selector.

Nevertheless, the manuscript’s central claim of a fully certified proof is not yet defensible. I found exact counterexamples to the stated outward-enclosure contracts, an incomplete rounding argument in the parameter selector, and verification claims that overstate independence and proof coverage. Additional hull-construction and resumption defects require repair or carefully proved restrictions.

**These findings do not disprove the conjecture or demonstrate that the final computation made an incorrect exclusion. They do invalidate parts of the justification offered for treating that computation as a proof.**

:chatgpt-content-reference{index="20"}[Download the reproducible referee audit, counterexamples and verification logs](sandbox:/mnt/data/mub6_referee_audit_2026-09-29.zip)

## 2. Scope and evidence limits

I reviewed the 17-page manuscript, its Markdown source, the computational modules, merged results, verification scripts, and accompanying logs, including the documented superseded runs. The standard applied is that of a research paper in mathematical quantum information or computational mathematics whose principal result depends on exhaustive, certified computation.

The audit distinguishes three kinds of evidence:

| Evidence class | What it establishes here |
|---|---|
| Mathematical argument | Whether the stated reductions and certificate implications follow from their hypotheses. |
| Executed checks | Whether specified data, algorithms and regression examples behave as reported. |
| Literature and editorial judgement | The contribution’s positioning, likely significance and prospective REF quality. |

I **did not rerun the complete phase-space branch-and-bound computation** reported as taking approximately 31 CPU-hours. I did complete the stored-certificate combinatorial check and a separate, vectorised bookkeeping check of the parameter cover.

The original Arb table verifier could not be executed because `python-flint` was unavailable. Instead, I checked **all 131,072 complex table entries—262,144 real components—using `mpmath` interval arithmetic at 40 decimal digits**. Every component satisfied the stated \(2^{-53}\) bound. Only after that check did I bypass the unavailable Arb dependency for the remaining tests. This is an independent software-level interval check, not formal verification of the interval library.

The archive was not text-indexed by the file-search service, so manuscript and code findings below are located by section, filename and function rather than invented file citations.

## 3. Contribution and positioning

The strongest defensible description is:

> The manuscript develops and executes a uniform parameter-box and phase-box exclusion scheme which, provided its numerical enclosure and coverage contracts are established, proves that no member of the closed Szöllősi family can belong to an MUB quadruplet.

That qualification is presently essential.

The mathematical organisation is good. The two-circulant family is represented through four phase variables subject to the complex equation
\[
x+y+\overline{xy}+u+v+\overline{uv}=0.
\]
The computation covers this constrained parameter set rather than relying on a generic choice of cubic roots. It then seeks enclosing phase boxes for all vectors unbiased to the standard basis and the selected Hadamard basis.

The combinatorial reduction is particularly clear. Given a covering collection of sufficiently small hulls, an MUB quadruplet would require twelve distinct hulls forming two six-cliques in the possible-orthogonality graph, with all cross-pairs present in the possible-unbiasedness graph. The exclusion is conservative: additional graph edges may make certification harder, but cannot invalidate the implication.

The intended result corresponds to Conjecture 3 in the published MMVW paper. Their Proposition 4.2 explains its relevance to the remaining algebraic Conjecture 2. The manuscript correctly retains that condition: **it does not claim an unconditional solution of the dimension-six MUB problem.** :chatgpt-content-reference{index="0"}

The likely originality lies in completing this particular full-family exclusion, together with its parameter-aware enclosure machinery and treatment of clustered solutions. The broad discretisation-and-combinatorics strategy has a substantial predecessor in Jaming et al.’s Fourier-family exclusion. It should not be presented as newly converting an otherwise merely numerical strategy into mathematics. :chatgpt-content-reference{index="1"}

## 4. Major comments

### 4.1. The generic outward-enclosure claim is false as implemented

**Issue.** Section 5.4 states that floating-point evaluation of the radii remains enclosing because the added error terms dominate the rounding error. That assertion is false for the supplied primitive operations.

**Manuscript evidence.** In `rig.py`, `dmul` computes the radius as
\[
|m_1|r_2+r_1|m_2|+r_1r_2+4U|m_1||m_2|,
\]
using ordinary floating-point arithmetic, without an outward-rounding step or a separate bound covering the evaluation of every radius term.

I executed the following exact regression. Set
\[
r=1+2^{-52}.
\]
The call
```python
rig.dmul(0j, r, 0j, r)
```
returns midpoint zero and radius
\[
\widehat r=1+2^{-51}.
\]
But the product of the input discs contains a point of modulus
\[
r^2=1+2^{-51}+2^{-104}.
\]
Thus the returned radius is too small by exactly
\[
2^{-104}\approx4.930380657631324\times10^{-32}.
\]

This is not a comparison against an approximate reference value: the shortfall was established with exact rational arithmetic.

A corresponding failure occurs in the one-by-one `imatmul` calculation. Importantly, it also occurs with **zero left radius**, taking
\[
A_m=r,\quad A_r=0,\qquad B_m=0,\quad B_r=r.
\]
Therefore, the zero-left-radius condition used at active call sites is not, by itself, enough to establish the general matrix-product contract.

**Why it matters.** The small size of the shortfall does not make a false enclosure statement true. The proof depends on enclosures being valid, not merely very accurate.

There is, however, an important limitation: **these examples do not establish that the final run encountered these inputs or that a downstream exclusion was wrong**. Additional slack or more restrictive reachable-input invariants might suffice. Those would need to be proved.

**Required revision.** Either implement rigorously outward evaluation of the radius arithmetic, or supply a complete error analysis for the actual reachable inputs and every final comparison. Audit the whole primitive layer, not just these two examples. Then rerun or independently replay every affected exclusion. A single final `nextafter` is not automatically sufficient for an arbitrarily long expression.

### 4.2. The parameter-selector exclusions lack the claimed certified arithmetic justification

**Issue.** The abstract and §5.1 say that no library trigonometric or exponential function enters the certifying path. The parameter selector contradicts that statement.

**Manuscript evidence.** `parambox2.py:test_boxes2` evaluates centre phases with `np.exp`. Its support-function calculation uses floating-point normalisation and square roots; the imported ordering and sector tests use floating-point angular constants. The module’s own justification is that errors are far below `EPS = 1e-9`.

The accompanying statement that candidate-direction rounding cannot affect soundness is incomplete. An arbitrarily chosen direction is acceptable, but its norm and the evaluated support-function inequality must still be controlled. A computed vector that is slightly longer than one, together with rounding in the inequality, requires an explicit bound.

**Why it matters.** This selector determines which parameter boxes are discarded before the phase-space certificates are built. An erroneous exclusion here would create a hole which no subsequent clique check could detect.

My cover audit confirms that the final results match **what this selector retains**. It does not prove that the selector retains every box meeting the required mathematical domain.

The large safety margin and conservative Taylor remainder make repair plausible. They are not substitutes for a complete argument.

**Required revision.** Provide either certified interval re-evaluation of all parameter exclusions or a documented end-to-end error bound covering centre evaluation, angular constants, normalisation, support-function arithmetic and comparison. Correct the abstract and §5.1 so that they accurately describe the implementation.

### 4.3. “Independent verification” is broader in the prose than in the supplied checks

**Issue.** The package contains valuable verification, but not an independent end-to-end validation of the entire proof.

**Manuscript evidence.**

`verify_cover.py` recomputes the cover using the original selector functions. This checks result-tree completeness and bookkeeping, but not the selector’s numerical soundness by an independent method.

`verify_cliques.py` uses a separate set-based clique enumeration, which is a real and useful independence improvement. However, it shares the supplied hull geometry, pair-interval construction and graph-building code.

The merged records store the final hulls rather than a complete trace of the phase-space subdivision, exclusion and contraction decisions. The **138 empty-hull records** contain no hull geometry from which an independent reader can reconstruct the phase-space exhaustion.

**Why it matters.** There are three separate questions:

1. Are all required parameter boxes represented?
2. Do the stored hulls cover every relevant solution?
3. Do those hulls exclude the required graph configuration?

The completed checks strongly support the first question relative to the selector, and the third relative to the numerical graph construction. They do not independently establish the second.

**Required revision.** Describe each verifier’s actual independence boundary. For end-to-end assurance, provide either an independently implemented certifier or a replayable subdivision/exclusion trace checked by a smaller outward-arithmetic kernel. A proof-assistant formalisation is not mandatory, but the claimed validation must extend beyond rechecking the final graphs.

### 4.4. `cyclic_hull` does not satisfy its documented general containment contract

**Issue.** The cyclic interval routine can return an interval that fails to contain its inputs.

**Manuscript evidence.** `pipeline.py:cyclic_hull` describes itself as producing a cyclic covering interval. On the period-\(65,536\) grid, give it
\[
[1000,1100],\qquad[65000,68000],
\]
repeated in the five coordinates. It returns
\[
[65000,66636].
\]
The input union contains residue \(2000\), whereas the returned hull does not.

**Why it matters.** Theorem 10 relies on the hull-construction step preserving all surviving solution boxes.

The qualification here is important. The active grouping path uses a narrower input pattern: a possibly wrapping declared ball together with non-wrapping retained dyadic boxes. The regression can omit part of the declared ball without proving that it omits an actual retained dyadic solution box. **I have not demonstrated loss of a solution on the executed path.**

Nevertheless, the documented generic contract is false, and the manuscript does not prove the narrower contract needed to bypass that failure.

**Required revision.** Repair the circular-union algorithm, for example by splitting wrapping intervals before constructing the complement gaps, or state and prove the exact restricted containment invariant used by the solver. Add adversarial wrapping and nested-interval tests; the supplied random tests did not expose this example.

### 4.5. Resumption can treat failed records as completed successes

**Issue.** The driver’s resumption semantics do not distinguish successful certification from a failed record awaiting subdivision.

**Manuscript evidence.** In a temporary results file, I inserted:
```json
{"lev": 9, "idx": [1, 2, 3, 4], "ok": false, "refined": true}
```
`driver.done_set` includes that key, and `driver.certify_box` returns one success and zero failures without verifying the descendants.

**Why it matters.** This is an operational correctness defect in a computation designed to support interruption, resumption and distributed execution.

It is **not evidence that the final merged run was contaminated**: every one of its 15,648 records has `ok=true`. The synthetic test exposes unsafe behaviour that matters particularly for the repair-and-rerun process.

**Required revision.** Separate “record exists”, “certified”, “failed”, and “subdivided with complete certified descendants”. Add crash-and-resume tests covering each state. A failed parent must never be counted as a successful certificate merely because its key exists.

### 4.6. Claims about prior rigour and exact degeneracy need narrowing

**Issue.** Two parts of the interpretation exceed their demonstrated support.

**Manuscript evidence.** Sections 1 and 10.2 contrast the present work with Jaming et al. in language suggesting that the earlier family exclusion was a numerical argument with an accuracy parameter. The earlier paper actually develops explicit mathematical discretisation bounds and a finite exclusion procedure. Its specific floating-point implementation may warrant scrutiny, but that requires a documented comparison rather than a blanket distinction between proof and numerics. :chatgpt-content-reference{index="2"}

Separately, §§8.2–8.3 use nearby numerical root clusters and small singular values to support statements that vectors coincide exactly and that solutions are not isolated.

**Why it matters.** A singular Jacobian does not, by itself, imply a non-isolated solution. Several branches may coalesce into an isolated multiple root. Nearby numerical clustering does not distinguish these cases.

Fortunately, the proposed covering certificate does not need the stronger claims.

**Required revision.** Present the degeneracy observations as numerical diagnostics unless an analytic argument establishes multiplicity or positive-dimensional solution sets. Replace the claim that “no isolation argument can succeed” with the narrower statement justified by the failure of the chosen simple-root isolation method. Position the computational improvements against the actual earlier proof strategy.

## 5. Rigour, results and inference

### Mathematical assessment

I found no counterexample to the normalisation argument, the Hadamard constraint, or the basic graph-certificate implication. Proposition 9 is a sound and useful reduction **conditional on its coverage and numerical-bound hypotheses**.

Several design choices are mathematically appropriate.

Numerical Newton searches need not find every solution, because their role is to group surviving boxes rather than decide which boxes survive. Likewise, an approximate inverse used as a Krawczyk preconditioner need not itself be a certified inverse: what matters is that the resulting operator is enclosed correctly. The proof also need not assume that the solution set is finite or consists only of simple roots.

The critical weakness is therefore not an obvious conceptual impossibility in the strategy. It is the gap between the stated numerical contracts and their implementation.

### Executed checks

| Check | Result | Interpretation |
|---|---:|---|
| Supplied checksum manifest | **63/63 matched** | Internal file consistency, not independent proof validity. |
| Complete unit table | **262,144 component checks; zero failures** | Supports the stated table-error bound. |
| Final records | **15,648 unique; all successful; all level 9** | Confirms the merged dataset’s stated structure. |
| Parameter bookkeeping | **422,400 descendants tested; 15,648 retained; zero missing or extra keys** | Complete relative to the supplied selector. |
| Stored graph certificates | **15,510 nonempty records; 6,628,140 cliques; zero reported problems** | Strong check of the combinatorial stage, with shared numerical geometry. |
| Numerical self-tests | **2,160 solution boxes; 34,560 pair checks; zero reported failures** | Useful consistency evidence, not exhaustive proof. |
| Symmetry tests | **200 random parameter points; deviations around \(10^{-15}\) or below** | Numerical support, not replacement for algebraic identities. |
| Adversarial regressions | **Enclosure, cyclic-hull and resumption defects reproduced** | Specific contracts require repair or narrower proof. |
| Complete phase-space proof search | **Not rerun** | No claim of independent end-to-end recertification. |

The successful reference clique enumeration is especially valuable given the previously discovered integer-bitmask bug. The manuscript deserves credit for reporting that failure, discarding the affected run and adding a differently implemented enumerator.

That history also illustrates why the new failures matter: ordinary successful runs and random self-tests do not establish universal software contracts.

## 6. External literature check

**Review date:** 29 September 2026.

I used targeted web searches, arXiv full-text records and version histories, the relevant journal publisher pages, and citation-chain checking. Representative search targets included “Szöllősi Conjecture 3 mutually unbiased bases”, the exact MMVW title, Jaming’s Fourier-family exclusion, the McNulty–Weigert comment and Liang–Chen reply, and recent order-six Hadamard classification claims. This was a targeted novelty and positioning check, not an exhaustive systematic review.

### Pivotal comparisons

**Szöllősi’s original family.** The original paper supplies the two-parameter construction and its hypocycloid-based parameter description. It supports the identification of the family under investigation; the present manuscript’s contribution is the proposed exclusion theorem, not discovery of that family. I consulted the full text. :chatgpt-content-reference{index="3"}

**Jaming et al.’s Fourier-family exclusion.** This is the main methodological antecedent. Its explicit discretisation inequalities and finite search show why the present work should identify its additional technical achievements precisely rather than claim ownership of the general strategy. I consulted the full text, including the discretisation argument. :chatgpt-content-reference{index="4"}

**MMVW, published in 2026.** The published article expressly formulates the target conjecture and anticipates a discretisation proof. This supports relevance and a genuine named target, but also means that implementing the proposed route is not, by itself, evidence of a wholly new conceptual method. I checked the published full text and the relevant conjecture/proposition. :chatgpt-content-reference{index="5"}

**The disputed earlier exclusion.** McNulty and Weigert identify a defect in an earlier lemma and its consequences for derived results. Liang and Chen’s latest arXiv version, dated **24 April 2025**, repairs one related statement but explicitly leaves the relevant lemma/theorem problem unresolved. The manuscript is therefore right not to regard the earlier family-exclusion proof as secure. Both texts were consulted; the repair should not be described as fixing nothing, nor as restoring the full exclusion. :chatgpt-content-reference{index="6"}

**The August 2026 classification claim.** Cárdenes Wuttig and Tindall’s order-six Hadamard classification preprint is already cited by the manuscript. I checked its abstract and version metadata, **not its full proof**. A Hadamard classification claim does not by itself verify this MUB exclusion, and should not be used as substitute evidence for it. :chatgpt-content-reference{index="7"}

Overall, the literature supports treating the target as meaningful. This search did not identify a verified duplicate establishing the same result by an already secure proof. That is a qualified novelty assessment, not a guarantee that no such work exists.

## 7. Minor comments

**Archive description.** The README says that three raw per-machine result files are omitted, while §9 describes them as supplied. Reconcile the statements. The merged file is sufficient for the checks I completed, but it is not identical to supplying the full raw provenance.

**Solution-count histogram.** The seven reported dominant categories total **15,580**, not 15,648. Another **68 records** have other numerical solution counts. Supply the complete histogram or label the listed categories as dominant rather than exhaustive. These remain numerical search counts, not certified exact solution counts.

**Figure interpretation.** A scatter plot of parameter-box centres is an illustration, not proof that the entire closed sector is covered. Attribute coverage to the verified selection/subdivision argument.

**Exact boundary example.** \(0.5+0.2887i\) is only approximately on the ray \(\arg\alpha=\pi/6\). An exact corresponding expression is
\[
\frac12+\frac{i}{2\sqrt3}.
\]

**Empty hulls.** The inference “every hull was excluded, i.e. the parameter box does not meet \(\Sigma\)” needs an additional existence argument. Empty common-unbiased-vector sets and absence of Hadamard parameters are not logically identical without that argument.

**Local consistency.** Reconcile descriptions of table-generation precision across documentation, and correct the module-count wording where “nine files” accompanies a ten-module table.

## 8. Prioritised revision plan

### Must fix before publishing the certified theorem

Establish sound outward arithmetic for every exclusion and contraction; certify the parameter selector’s complete error budget; repair or prove the restricted cyclic-hull invariant; and correct the resumption semantics. Then regenerate or independently replay the affected certificates. The existing final graphs alone cannot discharge these obligations.

Until this is done, describe the result as a **candidate computer-assisted proof**, not an established certified theorem.

### Should fix to strengthen the paper

Provide a precise dependency diagram showing what each verification script checks and what code it shares with the generator. Add replayable coverage/exclusion witnesses or an independent certifier. Include exact adversarial regression tests, reconcile the provenance inventory, and distinguish numerical degeneracy observations from analytic conclusions.

The methods comparison should identify exactly what improves on the earlier Fourier-family proof: parameter-aware contraction, clustered-solution handling, verification design, or measurable computational efficiency.

### Could improve presentation or future work

Extract a reusable theorem-and-checker interface for family-wide MUB exclusions. A clear specification of admissible families, enclosure obligations and graph witnesses would make the work easier to transfer beyond this one computation. This could strengthen its significance, but is not a prerequisite for publishing a correct full-family exclusion.

## 9. Editorial recommendation

**Major revisions.**

My confidence is **high** that the identified contract failures and reporting discrepancies are real. My confidence is **moderate** that they can be repaired without changing the headline conclusion: the successful checks and conservative design make that plausible, but plausibility is not certification.

I do not recommend rejection on grounds of a demonstrated false theorem. No such counterexample was found, and no irreparable mathematical obstruction to the proposed approach has been established.

Acceptance would require a sound arithmetic and coverage argument tied to an identified, successfully verified computation. Merely adjusting the prose, enlarging a few constants without proof, or rerunning the existing self-tests would not be sufficient.

## 10. Provisional REF calibration

The natural assessment home is **Unit of Assessment 10: Mathematical Sciences**, within Main Panel B. :chatgpt-content-reference{index="8"}

The current REF 2029 overview says that detailed panel criteria and working methods will be published in due course, including the quality-rating arrangements. I therefore use the established REF 2021 definitions as an **indicative calibration**, rather than claiming that final REF 2029 panel-specific criteria are available. In that framework, 4★ denotes world-leading quality, 3★ internationally excellent quality, and 2★ internationally recognised quality, assessed through originality, significance and rigour. :chatgpt-content-reference{index="9"}

| Dimension | As-submitted assessment | Rationale | Confidence |
|---|---|---|---|
| **Originality** | **3★** | A potentially new full-family exclusion with substantive implementation work; the overall strategy and suggested target are established antecedents. | Moderate |
| **Significance** | **3★** | Addresses a meaningful obstruction in the dimension-six programme, but does not establish the unconditional maximum-three result. | Moderate |
| **Rigour** | **2★** | Strong mathematical organisation and substantial reproducibility evidence, undermined by unresolved correctness contracts essential to the proof. | High on the defects; moderate on star calibration |
| **Overall output quality** | **2★** | The principal theorem cannot yet receive the credit warranted by an established proof. The rating is holistic, not an average of the rows. | Moderate |

### What successful revision could change

**A soundly repaired and reverified version would have a credible 3★ case.** It would establish a clearly identified theorem across a continuous family, including parameters that complicate simple-root methods.

**A 4★ case remains possible, but is not demonstrated by this submission.** A single family-exclusion theorem can be world-leading; there is no rule requiring a broader theorem. However, the case would need to rest on the result’s importance and the depth or transferability of the methods, supported by uncompromised verification—not simply on the number of boxes, the size of the computation, or the presence of a named conjecture.

My practical judgement for Evidence Press is therefore:

> **Retain this as a serious research candidate, prioritise the certification repairs, and do not currently present it as either a settled theorem or a four-star release.**

These grades are advisory judgements about output quality, not an official REF panel decision.

## 11. References

### Family, conjecture and computational antecedents

Jaming, P., Matolcsi, M., Móra, P., Szöllősi, F., & Weiner, M. (2009). A generalized Pauli problem and an infinite family of MUB-triplets in dimension 6. *Journal of Physics A: Mathematical and Theoretical, 42*(24), 245305. DOI: 10.1088/1751-8113/42/24/245305. :chatgpt-content-reference{index="10"}

Matolcsi, M., Matszangosz, Á. K., Varga, D., & Weiner, M. (2026). Triplets of mutually unbiased bases. *Journal of Algebraic Combinatorics, 63*, Article 26. DOI: 10.1007/s10801-026-01506-x. :chatgpt-content-reference{index="11"}

Szöllősi, F. (2010). A two-parameter family of complex Hadamard matrices of order 6 induced by hypocycloids. *Proceedings of the American Mathematical Society, 138*(3), 921–928. Full-text version: arXiv:0811.3930v2. :chatgpt-content-reference{index="12"}

### Corrections and recent developments

Cárdenes Wuttig, M., & Tindall, J. (2026). *A complete classification of complex Hadamard matrices of order six* (arXiv:2608.18053v2) [Preprint]. arXiv. DOI: 10.48550/arXiv.2608.18053. **Abstract and metadata checked only.** :chatgpt-content-reference{index="13"}

Liang, M., & Chen, L. (2025). *Reply to Comment on “Product states and Schmidt rank of mutually unbiased bases in dimension six”* (arXiv:2504.15576v3). arXiv. DOI: 10.48550/arXiv.2504.15576. :chatgpt-content-reference{index="14"}

McNulty, D., & Weigert, S. (2025). *Comment on “Product states and Schmidt rank of mutually unbiased bases in dimension six”* (arXiv:2504.13067). arXiv. DOI: 10.48550/arXiv.2504.13067. :chatgpt-content-reference{index="15"}

### Assessment guidance

Research Excellence Framework. (n.d.). *Guidance on REF 2021 results*. Official REF guidance, consulted 29 September 2026. :chatgpt-content-reference{index="16"}

Research Excellence Framework. (2025). *Section 1—Overview*. REF 2029 guidance, updated 10 December 2025; consulted 29 September 2026. :chatgpt-content-reference{index="17"}

Research Excellence Framework. (n.d.). *Main Panel B—Physical Sciences, Engineering and Mathematics*. REF 2029 panel information, consulted 29 September 2026. :chatgpt-content-reference{index="18"}