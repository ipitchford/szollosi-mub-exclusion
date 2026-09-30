"""Adversarial regression tests for the repaired GENERATOR (referee comments 4.1, 4.2, 4.4, 4.5).
Run from v2/:  python3 tests/test_gen.py"""
import sys, os, json, tempfile, numpy as np
from fractions import Fraction as Q
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gen'))
import rnd, core, hulls, selector, family, solve
fails = 0
def report(name, ok, detail=''):
    global fails
    print(f"[{'ok' if ok else 'FAIL'}] {name} {detail}", flush=True)
    if not ok: fails += 1

# ---- 4.1 referee counterexamples: disc product and point x interval matrix product
r = 1.0 + 2.0 ** -52
m, rad = core.dmul(np.array(0j), np.array(r), np.array(0j), np.array(r))
report("dmul encloses (referee r = 1+2^-52)", Q.from_float(float(rad)) >= Q.from_float(r) ** 2, f"radius {float(rad).hex()}")
Cm, Cr = core.imatmul_point(np.array([[r]]), np.array([[0.0]]), np.array([[r]]))
report("imatmul_point encloses (zero left radius)", Q.from_float(float(Cr[0, 0])) >= Q.from_float(r) ** 2, f"radius {float(Cr[0,0]).hex()}")
# randomised exact check of dmul against rational arithmetic (real inputs to keep the check exact)
rng = np.random.default_rng(0); bad = 0
for _ in range(20000):
    a, b = rng.uniform(-2, 2, 2); ra, rb = rng.uniform(0, 1, 2) * 10.0 ** rng.integers(-17, 1, 2)
    mm, rr = core.dmul(np.array(complex(a)), np.array(ra), np.array(complex(b)), np.array(rb))
    # worst case of |z1 z2 - mm| for z1 in [a-ra,a+ra], z2 in [b-rb,b+rb] (real): at the corners
    worst = max(abs(Q.from_float(x) * Q.from_float(y) - Q.from_float(float(mm.real)))
                for x in (a - ra, a + ra) for y in (b - rb, b + rb)) if False else None
    exact = max(abs((Q.from_float(a) + s * Q.from_float(ra)) * (Q.from_float(b) + t * Q.from_float(rb)) - Q.from_float(float(mm.real)))
                for s in (-1, 1) for t in (-1, 1))
    if Q.from_float(float(rr)) < exact: bad += 1
report("dmul encloses on 20,000 random real discs (exact rational check)", bad == 0, f"{bad} violations")

# ---- 4.4 cyclic hull: referee example and brute-force adversarial tests
N = 2 ** 16
lo = np.array([[1000] * 5, [65000] * 5]); hi = np.array([[1100] * 5, [68000] * 5])
hl, hh = hulls.cyclic_hull(lo, hi)
report("cyclic_hull referee example", bool(hulls.contains(hl, hh, lo, hi).all()), f"hull [{hl[0]},{hh[0]}]")
bad = 0; nonmin = 0
for trial in range(3000):
    k = int(rng.integers(1, 7)); M = [16, 64, 2 ** 16][trial % 3]
    a = rng.integers(-2 * M, 2 * M, k); ln = rng.integers(0, M, k)
    if trial % 5 == 0: ln[:] = rng.integers(0, 4, k)                     # tiny / nested
    if trial % 7 == 0: a[0] = M - 1 - ln[0] // 2                          # wrapping
    s, e = hulls.cyclic_hull_1d(a, a + ln, M)
    L = e - s
    for ai, li in zip(a, ln):
        off = (ai - s) % M
        if not (L >= M or off + li <= L): bad += 1
    # brute force: the hull may not be larger than the minimal covering arc on the integer grid (optimality)
    if M == 16:
        pts = set()
        for ai, li in zip(a, ln):
            for t in range(li + 1): pts.add((ai + t) % M)
        best = min(((st, ((max((p - st) % M for p in pts)))) for st in range(M)), key=lambda z: z[1]) if len(pts) < M else (0, M)
        if L > best[1] and not (L >= M): nonmin += 1
report("cyclic_hull contains all inputs (3,000 adversarial cases incl. wrapping/nested)", bad == 0, f"{bad} violations")
report("cyclic_hull minimal on small grids (informational)", nonmin == 0, f"{nonmin} non-minimal hulls")

# ---- 4.5 resumption states
import driver
with tempfile.TemporaryDirectory() as td:
    f = os.path.join(td, 'r.jsonl')
    open(f, 'w').write('\n'.join(json.dumps(x) for x in [
        {'lev': 9, 'idx': [1, 2, 3, 4], 'ok': False, 'status': 'refined'},
        {'lev': 9, 'idx': [5, 6, 7, 8], 'ok': False, 'status': 'failed'},
        {'lev': 9, 'idx': [1, 1, 1, 1], 'ok': True, 'status': 'certified'},
        {'lev': 9, 'idx': [2, 2, 2, 2], 'ok': True}]) + '\n')
    old = driver.RESULTS; driver.RESULTS = f
    try:
        s = driver.certified_set()
    finally:
        driver.RESULTS = old
report("resume: only 'certified' records are skipped", s == {(9, (1, 1, 1, 1))}, str(s))

# ---- 4.2 selector: every point of F lies in a kept box (sampled), and comparison with the old selector
from helpers_v2 import angles_for_alpha, in_region
miss = 0; tot = 0
while tot < 3000:
    a = complex(rng.uniform(0, 1.3), rng.uniform(0, 0.75))
    if not (0 <= np.angle(a) <= np.pi / 6) or not in_region(a): continue
    th = angles_for_alpha(a)
    for L in [7, 9, 11]:
        ii = np.floor(th / (2 * np.pi / 2 ** L)).astype(np.int64) % 2 ** L
        tot += 1
        if not selector.szollosi_test(ii[None, :], L)[0]: miss += 1
report("selector keeps boxes containing sampled points of F (levels 7, 9, 11)", miss == 0, f"{tot} boxes, {miss} missed")
# boundary cases: alpha on the ray and on the real axis, alpha = 0, alpha at the Dita point
for name, a in [('alpha=0', 0j), ('alpha=1', 1 + 0j), ('ray', 0.5 + 0.5j / np.sqrt(3)), ('real 0.7', 0.7 + 0j), ('ray 0.9', 0.9 * np.exp(1j * np.pi / 6))]:
    th = angles_for_alpha(a) if a != 0 else np.array([0, 2 * np.pi / 3, 0, 2 * np.pi / 3])
    ok = all(selector.szollosi_test((np.floor(th / (2 * np.pi / 2 ** L)).astype(np.int64) % 2 ** L)[None, :], L)[0] for L in [7, 9, 12])
    report(f"selector keeps the boxes of {name}", ok)
print(f"test_gen: {fails} failures"); sys.exit(1 if fails else 0)
