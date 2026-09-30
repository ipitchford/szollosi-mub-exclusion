"""probe_split_boxes.py -- NUMERICAL PROBE (not part of the proof): sample parameters on Sigma inside a certified box,
solve for the unbiased vectors numerically (Newton multistart) and count those outside the stored hulls.
Usage: python3 probe_split_boxes.py i1,i2,i3,i4 GEN.jsonl"""
import sys, json, numpy as np
sys.path.insert(0, 'gen')
import family, solve
fam = family.szollosi()
target = [int(t) for t in sys.argv[1].split(',')]
rec = None
for l in open(sys.argv[2]):
    r = json.loads(l)
    if r['lev'] == 9 and r['idx'] == target and r['status'] == 'certified': rec = r
hlo = np.array(rec['hull_lo']); hhi = np.array(rec['hull_hi'])
rng = np.random.default_rng(0)
lo = np.array(target) * 2 * np.pi / 512; hi = (np.array(target) + 1) * 2 * np.pi / 512
def inside(th):
    g = th / (2 * np.pi) * 65536
    for a, b in zip(hlo, hhi):
        L = b - a
        off = (g - a) % 65536
        if np.all((L >= 65536) | (off <= L)): return True
    return False
def surface_point(p):
    p = p.copy()
    for _ in range(60):
        x, y, u, v = np.exp(1j * p); z = np.conj(x * y); w = np.conj(u * v)
        C = x + y + z + u + v + w
        if abs(C) < 1e-14: break
        D = np.array([1j * (x - z), 1j * (y - z), 1j * (u - w), 1j * (v - w)])
        M = np.array([D.real, D.imag])                       # 2 x 4
        p = p - np.linalg.pinv(M) @ np.array([C.real, C.imag])   # minimum-norm projection onto Sigma
    return p, abs(C)
nout = 0; ntot = 0; worst = []
for trial in range(400):
    p0 = lo + rng.uniform(0, 1, 4) * (hi - lo)
    p, res = surface_point(p0)
    if res > 1e-12 or np.any(p < lo - 1e-12) or np.any(p > hi + 1e-12): continue
    sols = solve.multistart_solutions(fam.H(p), nstart=4000, seed=trial)
    x, y, u, v = np.exp(1j * p); tz = (-p[0] - p[1]) % (2 * np.pi); tw = (-p[2] - p[3]) % (2 * np.pi)
    inF = (p[0] <= p[1] <= tz) and (p[2] <= p[3] <= tw)
    out = [s for s in sols if not inside(s % (2 * np.pi))]
    ntot += len(sols); nout += len(out)
    if len(out) or trial % 40 == 0: print(f'p={np.round(p,5)} inF={inF} alpha={x+y+np.conj(x*y):.4f} sols={len(sols)} outside_hulls={len(out)}')
print('total', ntot, 'outside', nout)
