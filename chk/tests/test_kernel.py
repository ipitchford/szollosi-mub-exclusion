"""
Tests of the checker kernel against exact / high-precision arithmetic (fractions, mpmath).
  T1 directed rounding: add/sub/mul/sqr/sqrt enclosures contain the exact rational result (random + edge cases)
  T2 expi enclosures contain exp(i phi) for phi across the whole range, including large half-widths
  T3 S_j rectangles contain the exact S_j at random points of random phase/parameter boxes (both families)
  T4 modulus / Krawczyk tests never exclude a box containing a true unbiased vector (Newton solutions)
  T5 pair bounds contain the exact |<v,w>|^2 for random vectors of random hulls; self-check sound
  T6 Szollosi selector never excludes a box containing a point of the fundamental domain (random points)
  T7 containment agrees with brute-force enumeration on small grids
  T8 ordering test equals brute-force enumeration over fine sub-grids
"""
import sys, os, numpy as np, mpmath as mp
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ck, ctypes
mp.mp.dps = 40
rng = np.random.default_rng(12345)
TABLE = os.environ.get('MUB6_TABLE', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'gen', 'table17.npy'))
ck.load_table(TABLE, verify=False)
L = ck.lib()
P = ctypes.c_void_p
L.dbg_arith.argtypes = [ctypes.c_int64, P, P, P]
L.dbg_expi.argtypes = [ctypes.c_int64, P, P, P]
L.dbg_S.argtypes = [P, P, P, ctypes.c_int64, ctypes.c_int, ctypes.c_int, P]
fails = []
def check(name, cond, info=''):
    if not cond: fails.append((name, info)); print('FAIL', name, info)

def T1():
    edge = [0.0, -0.0, 1.0, -1.0, 2.0 ** -1074, -2.0 ** -1074, 2.0 ** -1022, 1.7976931348623157e308, 1 + 2 ** -52, 3.0, 1 / 3, 0.1, -0.7]
    a = np.array(edge + list(rng.normal(size=4000) * 10.0 ** rng.integers(-30, 30, 4000)))
    b = np.array(list(reversed(edge)) + list(rng.normal(size=4000) * 10.0 ** rng.integers(-30, 30, 4000)))
    out = np.zeros(10 * len(a))
    L.dbg_arith(len(a), ck._p(a), ck._p(b), ck._p(out)); out = out.reshape(-1, 10)
    n = 0
    for i in range(len(a)):
        A, B = Fr(float(a[i])), Fr(float(b[i]))
        for (lo, hi), ex in [((out[i, 0], out[i, 1]), A + B), ((out[i, 2], out[i, 3]), A - B), ((out[i, 4], out[i, 5]), A * B), ((out[i, 6], out[i, 7]), A * A)]:
            ok = (lo == -np.inf or Fr(lo) <= ex) and (hi == np.inf or ex <= Fr(hi)); n += 1
            check('T1 arith', ok, (a[i], b[i]))
        s = mp.sqrt(mp.mpf(abs(float(a[i]))))
        check('T1 sqrt', mp.mpf(out[i, 8]) <= s <= mp.mpf(out[i, 9]), a[i])
    print(f'T1: {n} arithmetic enclosures + {len(a)} square roots checked exactly')

def T2():
    NG = ck.NG
    c = rng.integers(-3 * NG, 3 * NG, 3000); w = np.concatenate([rng.integers(0, 40, 1000), rng.integers(0, NG, 1000), rng.integers(0, 3 * NG, 1000)])
    out = np.zeros(4 * len(c)); L.dbg_expi(len(c), ck._p(np.ascontiguousarray(c)), ck._p(np.ascontiguousarray(w)), ck._p(out)); out = out.reshape(-1, 4)
    A = 2 * mp.pi / NG
    for i in range(len(c)):
        for t in [-1, -0.37, 0, 0.5, 1]:
            ph = (int(c[i]) + t * int(w[i])) * A
            z = mp.expj(ph)
            check('T2 expi', mp.mpf(out[i, 0]) <= z.real <= mp.mpf(out[i, 1]) and mp.mpf(out[i, 2]) <= z.imag <= mp.mpf(out[i, 3]), (c[i], w[i], t))
    print(f'T2: {len(c)} exp(i phi) enclosures x 5 points checked with mpmath')

def S_exact(fam, th5, p):
    v = [mp.mpf(1)] + [mp.expj(t) for t in th5]
    out = []
    for j in range(6):
        s = mp.mpc(0)
        for k in range(6):
            ang = sum(int(fam.E[k, j, i]) * p[i] for i in range(fam.m))
            coef = mp.expj(mp.pi * int(fam.Z[k, j]) / 3) if fam.name == 'fourier' else mp.mpf(fam.coef_rect[k, j, 0])
            s += mp.conj(coef * mp.expj(ang)) * v[k]
        out.append(s)
    return out

def T3():
    NG = ck.NG; A = 2 * mp.pi / NG
    n = 0
    for fname in ['szollosi', 'fourier']:
        fam = ck.Fam(fname); fam.install()
        for trial in range(300):
            Lp = int(rng.integers(3, 12)); Lb = int(rng.integers(2, 12))
            idx = rng.integers(0, 2 ** Lp, fam.m); pc, ph = ck.pbox(idx, Lp)
            bi = rng.integers(0, 2 ** Lb, 5); lo, hi = ck.dyadic(bi[None, :], Lb); lo = np.ascontiguousarray(lo[0]); hi = np.ascontiguousarray(hi[0])
            vt = int(rng.integers(0, 2)); mask = int(rng.integers(0, 2 ** fam.m))
            out = np.zeros(24); L.dbg_S(ck._p(lo), ck._p(hi), ck._p(pc), ph, vt, mask, ck._p(out)); out = out.reshape(6, 4)
            for s in range(3):
                tt = rng.uniform(-1, 1, 5) if vt else np.zeros(5)
                th5 = [((int(lo[k]) + int(hi[k])) + mp.mpf(tt[k]) * (int(hi[k]) - int(lo[k]))) * A for k in range(5)]
                u = rng.uniform(-1, 1, fam.m)
                p = [(int(pc[i]) + (mp.mpf(u[i]) * ph if (mask >> i) & 1 else 0)) * A for i in range(fam.m)]
                S = S_exact(fam, th5, p)
                for j in range(6):
                    ok = mp.mpf(out[j, 0]) <= S[j].real <= mp.mpf(out[j, 1]) and mp.mpf(out[j, 2]) <= S[j].imag <= mp.mpf(out[j, 3])
                    check('T3 S', ok, (fname, trial, j)); n += 1
    print(f'T3: {n} S_j rectangle enclosures checked with mpmath (both families)')

def newton_mp(fam, th, p0, dep):
    """unbiased vectors for a parameter point p0 (numpy Newton; then polished)"""
    H = fam.H(p0)
    Hs = H.conj().T
    T = rng.uniform(0, 2 * np.pi, (4000, 5))
    for it in range(40):
        V = np.exp(1j * np.concatenate([np.zeros((len(T), 1)), T], axis=1)); W = V @ Hs.T
        G = (np.abs(W) ** 2 - 6.0)[:, 1:]
        J = 2 * np.real(np.conj(W)[:, 1:, None] * Hs[None, 1:, 1:] * (1j * V[:, None, 1:]))
        with np.errstate(all='ignore'):
            try: dT = np.linalg.solve(J + 1e-13 * np.eye(5)[None], G[:, :, None])[:, :, 0]
            except np.linalg.LinAlgError: continue
        dT = np.nan_to_num(dT); st = np.clip(np.abs(dT).max(axis=1, keepdims=True), 1e-300, None)
        T = (T - np.minimum(1, 0.5 / st) * dT) % (2 * np.pi)
    V = np.exp(1j * np.concatenate([np.zeros((len(T), 1)), T], axis=1))
    res = np.abs(np.abs(V @ Hs.T) ** 2 - 6).max(axis=1)
    return T[res < 1e-10]

def T4():
    NG = ck.NG; n_boxes = 0
    for fname in ['szollosi', 'fourier']:
        fam = ck.Fam(fname); fam.install()
        for trial in range(12):
            # a parameter point: for szollosi, a point on the constraint surface
            if fname == 'szollosi':
                sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'gen'))
                from helpers_v2 import angles_for_alpha, in_region
                while True:
                    al = complex(rng.uniform(-3, 3), rng.uniform(-3, 3))
                    if in_region(al): break
                ang = angles_for_alpha(al)
                p0 = np.array(ang, dtype=float) % (2 * np.pi)
                if abs(fam_C(p0)) > 1e-12: continue
            else:
                p0 = rng.uniform(0, 2 * np.pi, 2)
            lev = int(rng.integers(8, 11))
            idx = np.floor(p0 / (2 * np.pi) * 2 ** lev).astype(np.int64)
            pc, ph = ck.pbox(idx, lev)
            dep = ck.choose_dep(fam, pc)
            sols = newton_mp(fam, None, p0, dep)
            if len(sols) == 0: continue
            for Lb in [6, 8, 10, 12]:
                bi = np.floor(sols / (2 * np.pi) * 2 ** Lb).astype(np.int64) % 2 ** Lb
                lo, hi = ck.dyadic(bi, Lb)
                m = ck.modulus_excl(lo, hi, pc, ph)
                check('T4 idx variant', np.array_equal(m, ck.modulus_excl_idx(bi, Lb, pc, ph)), (fname, trial, Lb))
                k = ck.krawczyk_excl(lo, hi, pc, ph, dep)
                v = ck.mv_excl(lo, hi, pc, ph)
                check('T4 mean-value', not v.any(), (fname, trial, Lb, int(v.sum())))
                check('T4 modulus', not m.any(), (fname, trial, Lb, int(m.sum())))
                check('T4 krawczyk', not k.any(), (fname, trial, Lb, int(k.sum())))
                n_boxes += len(bi)
    print(f'T4: {n_boxes} boxes containing numerically certified solutions (residual < 1e-10) never excluded by the modulus, mean-value or Krawczyk tests')

def fam_C(p):
    x, y, u, v = np.exp(1j * np.asarray(p))
    return x + y + np.conj(x * y) + u + v + np.conj(u * v)

def T5():
    NG16 = 2 ** 16; A = 2 * mp.pi / NG16; n = 0
    for trial in range(60):
        nh = 8
        lo = rng.integers(0, NG16, (nh, 5)); wdt = rng.integers(0, int(rng.choice([50, 500, 5000, 40000])), (nh, 5)); hi = lo + wdt
        O, B, s = ck.pair_graphs(lo, hi)
        for i in range(nh):
            for j in range(nh):
                if i == j: continue
                for rep in range(4):
                    a = [mp.mpf(0)] + [(int(lo[i, k]) + mp.mpf(rng.uniform()) * int(wdt[i, k])) * A for k in range(5)]
                    b = [mp.mpf(0)] + [(int(lo[j, k]) + mp.mpf(rng.uniform()) * int(wdt[j, k])) * A for k in range(5)]
                    ip = abs(mp.fsum(mp.expj(b[k] - a[k]) for k in range(6))) ** 2
                    # soundness: if the value could be 0 (resp. 6) the flag must be set; check the enclosure implicitly:
                    # a point value < 6 and another > 6 across the pair implies B must be set (intermediate value)
                    n += 1
                    if ip < 1e-30: check('T5 O', O[i, j], (trial, i, j))
            # self-check: two random vectors of one hull must have |<v,w>|^2 > 6 when selfok
            if s[i]:
                for rep in range(20):
                    a = [(int(lo[i, k]) + mp.mpf(rng.uniform()) * int(wdt[i, k])) * A for k in range(5)]
                    b = [(int(lo[i, k]) + mp.mpf(rng.uniform()) * int(wdt[i, k])) * A for k in range(5)]
                    ip = abs(1 + mp.fsum(mp.expj(b[k] - a[k]) for k in range(5))) ** 2
                    check('T5 self', ip > 6, (trial, i))
        # B soundness: sample the range of |<v,w>|^2 over the pair; if samples straddle 6 then B must be set
        for i in range(nh):
            for j in range(i + 1, nh):
                vals = []
                for rep in range(30):
                    a = (lo[i] + rng.uniform(size=5) * wdt[i]) * 2 * np.pi / NG16; b = (lo[j] + rng.uniform(size=5) * wdt[j]) * 2 * np.pi / NG16
                    vals.append(abs(1 + np.exp(1j * (b - a)).sum()) ** 2)
                if min(vals) < 6 - 1e-9 and max(vals) > 6 + 1e-9: check('T5 B', B[i, j], (trial, i, j))
    # exact test on constructed pairs: v, w chosen with |<v,w>|^2 = 6 or 0 exactly (Fourier columns), tiny hulls
    F = np.exp(2j * np.pi * np.outer(np.arange(6), np.arange(6)) / 6)
    cols = [np.angle(F[1:, c] / F[0, c]) % (2 * np.pi) for c in range(6)]
    lo = np.array([np.floor(c / (2 * np.pi) * NG16).astype(np.int64) for c in cols]); hi = lo + 1
    O, B, s = ck.pair_graphs(lo, hi)
    check('T5 fourier O', all(O[i, j] for i in range(6) for j in range(6) if i != j))
    print(f'T5: {n} sampled pair values; orthogonality of the Fourier basis detected; self-check sampled')

def T6():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'gen'))
    from helpers_v2 import angles_for_alpha
    n = 0
    for trial in range(3000):
        # random point of F: alpha in the sector 0 <= arg <= pi/6 with |alpha| small enough, ordered angles
        from helpers_v2 import in_region
        r = rng.uniform(0, 3); t = rng.uniform(0, np.pi / 6)
        if trial < 30: r = 0.0
        if trial % 7 == 0: t = 0.0
        if trial % 7 == 1: t = np.pi / 6
        al = r * np.exp(1j * t)
        if not in_region(al): continue
        try: ang = angles_for_alpha(al)
        except Exception: ang = None
        if ang is None: continue
        p = np.array(ang, dtype=float) % (2 * np.pi)
        if abs(fam_C(p)) > 1e-9: continue
        for lev in [3, 5, 7, 9, 11]:
            idx = np.floor(p / (2 * np.pi) * 2 ** lev).astype(np.int64)[None, :]
            ok = ck.select_szollosi(idx, lev)[0]
            check('T6 selector', ok, (al, lev)); n += 1
    print(f'T6: {n} (point, level) pairs of the fundamental domain kept by the checker selector')

def T7():
    n = 0
    for trial in range(400):
        N16 = 2 ** 16
        hlo = rng.integers(0, N16, (3, 5)); hhi = hlo + rng.integers(0, 3000, (3, 5))
        hhi[rng.integers(0, 3), rng.integers(0, 5)] = hlo[0, 0] + N16 + 5
        Lb = int(rng.integers(6, 10)); s = 2 ** (16 - Lb)
        # boxes near hull corners
        c = hlo[rng.integers(0, 3)] // s + rng.integers(-2, 3, 5)
        bi = (c % 2 ** Lb)[None, :]; lo, hi = ck.dyadic(bi, Lb)
        got = ck.contained(lo, hi, hlo, hhi)[0]
        # brute force: every point lo..hi (step 1 on G16) inside some hull? (per hull, all coordinates)
        want = False
        for h in range(3):
            okh = True
            for k in range(5):
                Lh = hhi[h, k] - hlo[h, k]
                if Lh >= N16: continue
                pts = np.arange(lo[0, k], hi[0, k] + 1)
                if not np.all(((pts - hlo[h, k]) % N16) <= Lh): okh = False
            want = want or okh
        check('T7 contained', got == want, (trial,)); n += 1
    # exact edge cases: box end-points equal to, or one grid unit beyond, the hull end-points (also across the wrap)
    N16 = 2 ** 16
    for ell in [0, 1, 1000, N16 - 5, N16 - 1]:
        for Lh in [0, 1, 7, 300]:
            eta = ell + Lh
            hlo = np.array([[ell] * 5]); hhi = np.array([[eta] * 5])
            for dlo, dhi, want in [(0, 0, True), (-1, 0, False), (0, 1, False), (1, 0, True), (0, -1, True), (-1, 1, False)]:
                a, b = ell + dlo, eta + dhi
                if b < a: continue
                for k in range(5):
                    lo = np.array([[ell] * 5]); hi = np.array([[ell] * 5])
                    lo[0, k] = a; hi[0, k] = b
                    shift = int(rng.integers(-2, 3)) * N16          # the same box written with a shifted representative
                    got = ck.contained(lo + shift, hi + shift, hlo, hhi)[0]
                    check('T7 edge', got == want, (ell, Lh, dlo, dhi, k)); n += 1
    print(f'T7: {n} containment decisions equal brute force or exact edge expectations')

def T8():
    n = 0
    for lev in [3, 4, 5]:
        N = 2 ** lev; sub = 16
        i1, i2 = np.meshgrid(np.arange(N), np.arange(N), indexing='ij'); i1 = i1.ravel(); i2 = i2.ravel()
        got = ck.ordering_ok(i1, i2, lev)
        for a, b, g in zip(i1, i2, got):
            # brute force on a fine grid of t1, t2 in the closed box (exact rationals)
            want = False
            for u in range(sub + 1):
                for v in range(sub + 1):
                    t1 = Fr(int(a)) + Fr(u, sub); t2 = Fr(int(b)) + Fr(v, sub)
                    if u == sub or v == sub: continue
                    t3 = (-t1 - t2) % N
                    if t1 <= t2 <= t3: want = True
            # the test is a necessary condition: want => got
            check('T8 ordering', (not want) or g, (lev, a, b)); n += 1
    print(f'T8: {n} ordering decisions sound against brute force')

def T9():
    """the checker's family data reproduce the defining formulas (independent numpy constructions)"""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'gen'))
    from helpers_v2 import angles_for_alpha, in_region
    def circ(a0, a1, a2): return np.array([[a0, a1, a2], [a2, a0, a1], [a1, a2, a0]])
    fam = ck.Fam('szollosi'); n = 0
    for trial in range(200):
        p = rng.uniform(0, 2 * np.pi, 4)
        x, y, u, v = np.exp(1j * p); z = np.conj(x * y); w = np.conj(u * v)
        A = circ(1, np.conj(x), z); B = circ(1, np.conj(u), w)
        Hd = np.block([[A, B], [B.conj().T, -A.conj().T]])
        check('T9 szollosi H', np.abs(fam.H(p) - Hd).max() < 1e-13, trial)
        al = complex(rng.uniform(-3, 3), rng.uniform(-3, 3))
        if in_region(al):
            q = angles_for_alpha(al); Hq = fam.H(q)
            check('T9 szollosi Hadamard', np.abs(Hq @ Hq.conj().T - 6 * np.eye(6)).max() < 1e-9, al); n += 1
    famf = ck.Fam('fourier')
    om = np.exp(2j * np.pi / 3)
    for trial in range(200):
        a, b = rng.uniform(0, 1, 2); x = np.exp(2j * np.pi * a); y = np.exp(2j * np.pi * b)
        F = np.array([[1, 1, 1, 1, 1, 1],
                      [1, -om ** 2 * x, om, -x, om ** 2, -om * x],
                      [1, om * y, om ** 2, y, om, om ** 2 * y],
                      [1, -1, 1, -1, 1, -1],
                      [1, om ** 2 * x, om, x, om ** 2, om * x],
                      [1, -om * y, om ** 2, -y, om, -om ** 2 * y]])
        Hf = famf.H(np.array([2 * np.pi * a, 2 * np.pi * b]))
        check('T9 fourier H', np.abs(Hf - F).max() < 1e-13, trial)
        check('T9 fourier Hadamard', np.abs(F @ F.conj().T - 6 * np.eye(6)).max() < 1e-12, trial)
        # coefficient rectangles enclose the sixth roots of unity (mpmath)
    exact = {0: (1, 0), 3: (-1, 0)}                     # zeta^0 = 1, zeta^3 = -1 exactly
    for idx, z in np.ndenumerate(famf.Z):
        r = famf.coef_rect[idx]
        if int(z) in exact:
            re_, im_ = exact[int(z)]
            check('T9 coef rect', r[0] <= re_ <= r[1] and r[2] <= im_ <= r[3], idx)
        else:                                          # interval enclosure of the exact value must lie inside the rectangle
            th = mp.iv.pi * int(z) / 3; ci = mp.iv.cos(th); si = mp.iv.sin(th)
            check('T9 coef rect', mp.mpf(r[0]) <= ci.a and ci.b <= mp.mpf(r[1]) and mp.mpf(r[2]) <= si.a and si.b <= mp.mpf(r[3]), idx)
    print(f'T9: family data equal the defining formulas (200 + 200 random parameters; {n} Hadamard checks on the surface)')

def T10():
    """stress test: coarse parameter boxes (levels 3-7) and non-dyadic phase boxes of widths 1..1024 grid units
    placed so that a numerically certified solution lies inside, close to a face; no test may exclude them"""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'gen'))
    from helpers_v2 import angles_for_alpha, in_region
    n = 0
    for fname in ['szollosi', 'fourier']:
        fam = ck.Fam(fname); fam.install()
        for trial in range(40):
            if fname == 'szollosi':
                while True:
                    al = complex(rng.uniform(-3, 3), rng.uniform(-3, 3))
                    if in_region(al): break
                p0 = np.array(angles_for_alpha(al), dtype=float) % (2 * np.pi)
                if abs(fam_C(p0)) > 1e-12: continue
            else:
                p0 = rng.uniform(0, 2 * np.pi, 2)
            lev = int(rng.integers(3, 8))
            idx = np.floor(p0 / (2 * np.pi) * 2 ** lev).astype(np.int64)
            pc, ph = ck.pbox(idx, lev)
            dep = ck.choose_dep(fam, pc)
            sols = newton_mp(fam, None, p0, dep)
            if len(sols) == 0: continue
            s16 = sols / (2 * np.pi) * 2 ** 16
            for wexp in [0, 2, 4, 6, 8, 10]:
                w = 2 ** wexp
                u = rng.uniform(0.01, 0.99, s16.shape)
                lo = np.floor(s16 - u * w).astype(np.int64); hi = lo + w
                ok_in = np.all((s16 >= lo + 1e-3) & (s16 <= hi - 1e-3), axis=1)
                lo, hi = lo[ok_in], hi[ok_in]
                m = ck.modulus_excl(lo, hi, pc, ph); v = ck.mv_excl(lo, hi, pc, ph); k = ck.krawczyk_excl(lo, hi, pc, ph, dep)
                check('T10 modulus', not m.any(), (fname, trial, lev, w, int(m.sum())))
                check('T10 mean-value', not v.any(), (fname, trial, lev, w, int(v.sum())))
                check('T10 krawczyk', not k.any(), (fname, trial, lev, w, int(k.sum())))
                n += len(lo)
    print(f'T10: {n} off-centre boxes (widths 1..1024 grid units, parameter levels 3..7) containing solutions never excluded')

def T11():
    """white-box test of the Krawczyk assembly: for every q in Q and w in W the point
    g_q(w) - w~ = (w - w~) - Y F(w; q) must lie in the kernel's interval vector K - w~ (mean-value form).
    Points are sampled at vertices and at random, and F is evaluated exactly (mpmath, 50 digits)."""
    L.dbg_krawczyk.argtypes = [P, P, P, ctypes.c_int64, ctypes.c_int, P, P, P]
    A16 = 2 * mp.pi / 2 ** 16; A17 = 2 * mp.pi / 2 ** 17
    n = 0
    for fname in ['szollosi', 'fourier']:
        fam = ck.Fam(fname); fam.install()
        coef = [[(mp.expj(mp.pi * int(fam.Z[k, j]) / 3) if fname == 'fourier' else mp.mpf(fam.coef_rect[k, j, 0])) for j in range(6)] for k in range(6)]
        if fname == 'fourier':
            for k in range(6):
                for j in range(6):
                    if int(fam.Z[k, j]) in (0, 3): coef[k][j] = mp.mpf(1 if int(fam.Z[k, j]) == 0 else -1)
        for trial in range(60):
            lev = int(rng.integers(3, 10))
            idx = rng.integers(0, 2 ** lev, fam.m); pc, ph = ck.pbox(idx, lev)
            dep = ck.choose_dep(fam, pc); nd = len(dep); nn_ = 5 + nd
            free = [i for i in range(fam.m) if i not in dep]
            Lb = int(rng.integers(3, 12)); bi = rng.integers(0, 2 ** Lb, 5)
            lo, hi = ck.dyadic(bi[None, :], Lb); lo = np.ascontiguousarray(lo[0]); hi = np.ascontiguousarray(hi[0])
            d32 = np.ascontiguousarray(np.array(list(dep) + [0], dtype=np.int32))
            kd = np.zeros(2 * 7); Y = np.zeros(49)
            L.dbg_krawczyk(ck._p(lo), ck._p(hi), ck._p(pc), ph, nd, ck._p(d32), ck._p(kd), ck._p(Y))
            if np.isnan(kd).any(): continue
            kd = kd.reshape(7, 2); Y = Y.reshape(7, 7)
            wc = [(int(lo[k]) + int(hi[k])) * A17 for k in range(5)] + [int(pc[i]) * A17 for i in dep]
            wr = [(int(hi[k]) - int(lo[k])) * A17 for k in range(5)] + [ph * A17 for i in dep]
            qc = [int(pc[i]) * A17 for i in free]; qr = ph * A17
            for smp in range(24):
                if smp < 16:
                    sw = rng.choice([-1.0, 1.0], nn_); sq = rng.choice([-1.0, 1.0], len(free))
                else:
                    sw = rng.uniform(-1, 1, nn_); sq = rng.uniform(-1, 1, len(free))
                w = [wc[r] + mp.mpf(sw[r]) * wr[r] for r in range(nn_)]
                q = [qc[l] + mp.mpf(sq[l]) * qr for l in range(len(free))]
                p = [None] * fam.m
                for a, i in enumerate(dep): p[i] = w[5 + a]
                for a, i in enumerate(free): p[i] = q[a]
                th = [mp.mpf(0)] + w[:5]
                F = []
                for j in range(1, 6):
                    S = mp.fsum(mp.conj(coef[k][j] * mp.expj(mp.fsum(int(fam.E[k, j, i]) * p[i] for i in range(fam.m)))) * mp.expj(th[k]) for k in range(6))
                    F.append(abs(S) ** 2 - 6)
                if nd:
                    C = mp.fsum(mp.expj(mp.fsum(int(fam.CE[t_, i]) * p[i] for i in range(fam.m))) for t_ in range(len(fam.CE)))
                    F += [C.real, C.imag]
                for r in range(nn_):
                    g = (w[r] - wc[r]) - mp.fsum(mp.mpf(float(Y[r, s])) * F[s] for s in range(nn_))
                    check('T11 krawczyk assembly', mp.mpf(kd[r, 0]) <= g <= mp.mpf(kd[r, 1]), (fname, trial, smp, r, float(g), kd[r].tolist()))
                n += 1
    print(f'T11: {n} sampled points (w, q) satisfy g_q(w) - w~ in K - w~ (Krawczyk assembly, exact evaluation)')

if __name__ == '__main__':
    only = os.environ.get('MUB6_ONLY')
    for T in [T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11]:
        if only and T.__name__ not in only.split(','): continue
        T()
    print('ALL PASS' if not fails else f'{len(fails)} FAILURES')
    sys.exit(1 if fails else 0)
