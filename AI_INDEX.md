# AI index: Szöllősi-family MUB exclusion

## Entry points

- Statement and proof: [paper/paper.md](paper/paper.md),
  [PDF](paper/paper.pdf), [TeX](paper/paper.tex).
- Current status: [STATUS.md](STATUS.md), [local replay receipt](REPLAY_RECEIPT.md).
- Trust boundary: [ASSURANCE.md](ASSURANCE.md).
- Reproduction contract: [CERTIFICATE_INTERFACE.md](CERTIFICATE_INTERFACE.md).
- Provenance: [PROVENANCE.md](PROVENANCE.md).
- Source bibliography: [paper/sections/partS.md](paper/sections/partS.md).

## Exact scope

For every Hadamard matrix in the specified closed two-circulant Szöllősi family,
the standard basis and its normalized columns cannot be extended by two further
mutually unbiased bases in complex dimension six. Boundary and repeated-root
parameters are included. The Fourier-family exclusion is a re-proof of prior
work. The unrestricted MUB(6) maximum remains conditional on an unproved
external conjecture; Karlsson's larger family is not covered.

## Evidence map

1. `paper/sections/part2.md`: abstract certificate theorem and obligations.
2. `paper/sections/part3.md`: family, symmetry reduction, owner-retaining selector.
3. `chk/ck.py`, `chk/kernel.c`: interval and graph checker.
4. `chk/run_check.py`, `chk/verify_split.py`: cover and subdivision drivers.
5. `chk/receipt_binding.py`: exact input/environment binding, not authentication.
6. `chk/full_replay.py`: clean two-family replay without old success receipts.
7. `chk/tests/test_revision.py`: regression and semantic negative controls.
8. `results/`: untrusted hull proposals; `gen/table17.npy`: verified table input.
9. `chk_results/`: HISTORICAL evidence only; consult current release receipt
   before treating any result as a replay of the revised source.
10. [Szöllősi complete summary](replay/szollosi_summary.json) and
    [Fourier complete summary](replay/fourier_summary.json): fresh local replay,
    exact source/data/binary/runtime identities. These are not authenticated
    external verification or a substitute for executing the clean driver.
11. [Claims](CLAIMS.json), [prior-art and limitation record](RESEARCH_GATES.json),
    [licence map](LICENSES.md) and [citation metadata](CITATION.cff).

## Safe reuse

Read the theorem hypotheses, then execute the clean command in the interface.
Require both complete covers and preserve environment-bound outputs. Graph
checks alone do not prove phase-space coverage. Hull counts are not root counts.
Do not infer external reproduction, formal verification, novelty or the full
dimension-six conjecture from a successful build, DOI or archive manifest.
