"""inventory.py -- Appendix C: inventory of code and data files (lines, SHA-256) for the paper."""
import os, hashlib, sys, json
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
groups = [
    ('Checker (trusted)', ['chk/kernel.c', 'chk/ck.py', 'chk/run_check.py', 'chk/verify_split.py', 'chk/receipt_binding.py', 'chk/full_replay.py', 'chk/build.sh']),
    ('Exact verification (trusted)', ['exact/lemma5_exact.py', 'fourier/symmetries.py', 'fourier/fundamental_domain.py', 'fourier/fourier_exact.py']),
    ('Tests', ['chk/tests/test_kernel.py', 'chk/tests/test_revision.py', 'chk/tests/test_boundary_interval.py', 'chk/tests/mutate.py', 'tests/test_gen.py']),
    ('Generator (untrusted)', ['gen/rnd.py', 'gen/family.py', 'gen/core.py', 'gen/hulls.py', 'gen/comb.py', 'gen/selector.py', 'gen/solve.py', 'gen/driver.py', 'gen/helpers_v2.py']),
    ('Diagnostics and analysis (not part of the proof)', ['exact/degeneracy_diagnostics.py', 'analysis/stats.py', 'analysis/empty_boxes.py', 'analysis/ready_roots.py', 'analysis/probe_split_boxes.py', 'analysis/special_boxes.py']),
    ('Data', ['gen/table17.npy']),
]
def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()
lines = ['# Appendix C. Inventory of code and data', '',
         'Line counts and SHA-256 digests (first 16 hexadecimal digits) of the files on which the results depend; the bundle file `SHA256SUMS` lists the full digests of every file, including all results, checker outputs and logs.', '']
tot = {}
for g, fs in groups:
    lines += [f'**{g}**', '', '| file | lines | SHA-256 (prefix) |', '|---|---:|---|']
    for f in fs:
        p = os.path.join(ROOT, f)
        n = sum(1 for _ in open(p, 'rb')) if not f.endswith('.npy') else '—'
        lines.append(f'| `{f}` | {n} | `{sha(p)[:16]}` |')
        if isinstance(n, int): tot[g] = tot.get(g, 0) + n
    lines.append('')
lines += ['Totals: ' + '; '.join(f'{g.lower()}, {n} lines' for g, n in tot.items()) + '.', '',
          '*Versions actually used.* Every checker run used `kernel.c`, `ck.py` and `run_check.py` with the digests above (they are recorded in each checker summary). The parameter-subdivision verifications used in the final cover evaluation were all run with the `verify_split.py` listed above; three earlier subdivision runs on machine B, with a version that differed only in one logging statement, are kept in `chk_results/superseded_split_runs/` and are not used. The generator code was identical on both machines except `driver.py` (machine A\'s copy has the tail-first option) and `hulls.py`, whose SciPy call was replaced during the run (Section 9.4); since the generator is not trusted, these differences do not affect the proof.', '']
lines[-2] = '*Version boundary.* These hashes describe revised source. Historical logs refer to earlier versions and are not fresh receipts for these bytes. Revised clean-replay summaries bind their exact sources, compiled kernel, certificate inputs and environment. Generator variations cannot replace checker verification.'
open(os.path.join(ROOT, 'paper', 'sections', 'partC.md'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
