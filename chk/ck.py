"""
ck.py -- INDEPENDENT CHECKER (trusted component; shares no code with ../gen).

Given a family, the list of certified parameter boxes and, for each box, the hulls produced by the (untrusted)
generator, the checker establishes with its own code and arithmetic:
  (1) TABLE    : every entry of the unit table encloses exp(2 pi i m / 2^17) within 2^-53 per component
                 (mpmath interval arithmetic -- a different library from the generator's Arb check);
  (2) COVER    : every point of the fundamental domain lies in a box selected by the checker's own selector,
                 and every selected box is certified (or proved, by deeper subdivision, to miss the domain);
  (3) COVERAGE : for every certified box P, every vector unbiased to I and H(p) (p in P, constraint holding)
                 lies in one of the box's hulls -- by an exhaustive branch-and-bound over the phase torus in
                 which a phase box is discarded only if the C kernel's modulus test or parametric Krawczyk test
                 excludes it, or if it is contained in a hull;
  (4) GRAPHS   : pair bounds of the hulls (C kernel), the self-check, and the absence of two 6-cliques of the
                 possibly-orthogonal graph joined across the possibly-unbiased graph (own clique search).
The C kernel (kernel.c) does all floating-point work with one-ulp outward rounding of every operation.
"""
import ctypes, os, sys, json, time, numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NG = 2 ** 17; NG16 = 2 ** 16
_lib = None

def lib():
    global _lib
    if _lib is None:
        alt = os.environ.get('MUB6_KERNEL')          # test-suite only (mutation testing); never set in a proof run
        names = ('libmubk.dylib',) if sys.platform == 'darwin' else ('libmubk.so',)
        for p in ([alt] if alt else [os.path.join(HERE, n) for n in names]):
            if os.path.exists(p):
                _lib = ctypes.CDLL(p); break
        if _lib is None: raise RuntimeError('build the kernel first (build.sh)')
        P = ctypes.c_void_p; I64 = ctypes.c_int64; I = ctypes.c_int; D = ctypes.c_double
        _lib.set_table.argtypes = [P, D]
        _lib.set_family.argtypes = [I, P, P, I, P, P]
        _lib.modulus_test.argtypes = [I64, P, P, P, I64, P]
        _lib.modulus_test_idx.argtypes = [I64, P, I, P, I64, P]
        _lib.krawczyk_test.argtypes = [I64, P, P, P, I64, I, P, P]
        _lib.mv_test.argtypes = [I64, P, P, P, I64, P]
        _lib.contained.argtypes = [I64, P, P, I64, P, P, P]
        _lib.pair_bounds.argtypes = [I64, P, P, P, P, P]
        _lib.szollosi_param_test.argtypes = [I64, P, I, D, D, P, P]
    return _lib

def _p(a): return a.ctypes.data_as(ctypes.c_void_p)

# ------------------------------------------------------------------------------------------------ (1) table
def verify_table(T, log=print):
    """every component within 2^-53 of the exact value, checked with mpmath.iv (40 digits)"""
    import mpmath as mp
    mp.iv.dps = 40
    bound = mp.iv.mpf(2) ** -53
    bad = 0; t0 = time.time()
    for m in range(len(T)):
        th = 2 * mp.iv.pi * m / len(T)
        if not bool(abs(mp.iv.cos(th) - mp.iv.mpf(float(T[m].real))).b <= bound): bad += 1
        if not bool(abs(mp.iv.sin(th) - mp.iv.mpf(float(T[m].imag))).b <= bound): bad += 1
    log(f"table: {len(T)} entries, {2*len(T)} components checked with mpmath.iv, {bad} failures [{time.time()-t0:.0f}s]")
    return bad == 0

def load_table(path, verify=True, log=print):
    T = np.load(path)
    if T.shape != (NG,) or T.dtype != np.complex128:
        raise ValueError('invalid unit table shape or dtype')
    if verify and not verify_table(T, log): raise RuntimeError('table verification failed')
    flat = np.ascontiguousarray(np.stack([T.real, T.imag], axis=1).reshape(-1))
    pi_up = float(np.nextafter(np.pi, np.inf))
    import mpmath as mp
    with mp.workdps(60):
        if not mp.mpf(np.pi) < +mp.pi < mp.mpf(pi_up):
            raise RuntimeError('pi enclosure failed')
    lib().set_table(_p(flat), pi_up)
    return T

# ------------------------------------------------------------------------------------------------ families
class Fam:
    """checker's own family description, built from the defining formulas"""
    def __init__(self, name):
        self.name = name
        if name == 'szollosi':
            self.m = 4
            # A = circ(1, conj x, conj(xy)), B = circ(1, conj u, conj(uv)); H = [A B; B* -A*]
            # exponent vectors of the entries a..f as functions of (th_x, th_y, th_u, th_v)
            ea = [0, 0, 0, 0]; eb = [-1, 0, 0, 0]; ec = [-1, -1, 0, 0]
            ed = [0, 0, 0, 0]; ee = [0, 0, -1, 0]; ef = [0, 0, -1, -1]
            neg = lambda v: [-t for t in v]
            def circ(p, q, r): return [[p, q, r], [r, p, q], [q, r, p]]
            A = circ(ea, eb, ec); Bm = circ(ed, ee, ef)
            Bs = [[neg(Bm[j][i]) for j in range(3)] for i in range(3)]      # B* : conjugate transpose
            As = [[neg(A[j][i]) for j in range(3)] for i in range(3)]
            E = np.zeros((6, 6, 4), dtype=np.int32); sign = np.ones((6, 6))
            for i in range(3):
                for j in range(3):
                    E[i, j] = A[i][j]; E[i, j + 3] = Bm[i][j]; E[i + 3, j] = Bs[i][j]; E[i + 3, j + 3] = As[i][j]
                    sign[i + 3, j + 3] = -1.0
            self.E = E
            self.coef_rect = np.zeros((6, 6, 4)); self.coef_rect[..., 0] = sign; self.coef_rect[..., 1] = sign
            self.CE = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [-1, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, -1, -1]], dtype=np.int32)
            self.ct_rect = np.tile(np.array([1.0, 1.0, 0.0, 0.0]), (6, 1))
        elif name == 'fourier':
            self.m = 2
            # sqrt6 F(a,b), Jaming et al. (2009) eq. (3); w = exp(2 pi i/3); x = e^{is}, y = e^{it}
            # entry = zeta^Z * x^X * y^Y, zeta = exp(i pi/3)
            rows = [["1", "1", "1", "1", "1", "1"],
                    ["1", "-w2 x", "w", "-x", "w2", "-w x"],
                    ["1", "w y", "w2", "y", "w", "w2 y"],
                    ["1", "-1", "1", "-1", "1", "-1"],
                    ["1", "w2 x", "w", "x", "w2", "w x"],
                    ["1", "-w y", "w2", "-y", "w", "-w2 y"]]
            Z = np.zeros((6, 6), dtype=int); E = np.zeros((6, 6, 2), dtype=np.int32)
            for k in range(6):
                for j in range(6):
                    s = rows[k][j]; z = 0
                    if s.startswith('-'): z += 3; s = s[1:]
                    toks = s.split()
                    for tk in toks:
                        if tk == 'w': z += 2
                        elif tk == 'w2': z += 4
                        elif tk == 'x': E[k, j, 0] = 1
                        elif tk == 'y': E[k, j, 1] = 1
                    Z[k, j] = z % 6
            self.E = E; self.Z = Z
            self.coef_rect = sixth_root_rects(Z)
            self.CE = np.zeros((0, 2), dtype=np.int32); self.ct_rect = np.zeros((0, 4))
        else:
            raise ValueError(name)
        self.nt = len(self.CE)
    def install(self):
        E = np.ascontiguousarray(self.E.reshape(-1).astype(np.int32))
        cr = np.ascontiguousarray(self.coef_rect.reshape(-1).astype(float))
        CE = np.ascontiguousarray(self.CE.reshape(-1).astype(np.int32))
        ct = np.ascontiguousarray(self.ct_rect.reshape(-1).astype(float))
        self._keep = (E, cr, CE, ct)
        if lib().set_family(self.m, _p(E), _p(cr), self.nt, _p(CE), _p(ct)) != 0:
            raise RuntimeError('kernel rejected family')
    def H(self, p):
        """numerical matrix (tests only)"""
        coef = np.exp(1j * np.pi * self.Z / 3) if self.name == 'fourier' else self.coef_rect[..., 0]
        return coef * np.exp(1j * (self.E @ np.asarray(p, dtype=float)))

def sixth_root_rects(Z):
    """rectangles [re_lo, re_hi, im_lo, im_hi] enclosing zeta^Z exactly (mpmath interval arithmetic)"""
    import mpmath as mp
    mp.iv.dps = 30
    out = np.zeros(Z.shape + (4,))
    for idx, z in np.ndenumerate(Z):
        th = mp.iv.pi * int(z) / 3
        c = mp.iv.cos(th); s = mp.iv.sin(th)
        out[idx] = [float(np.nextafter(float(c.a), -np.inf)), float(np.nextafter(float(c.b), np.inf)),
                    float(np.nextafter(float(s.a), -np.inf)), float(np.nextafter(float(s.b), np.inf))]
        if z == 0: out[idx] = [1.0, 1.0, 0.0, 0.0]
        if z == 3: out[idx] = [-1.0, -1.0, 0.0, 0.0]
    return out

# ------------------------------------------------------------------------------------------------ kernel wrappers
def pbox(idx, lev):
    ph = 2 ** (16 - lev)
    pc = np.ascontiguousarray(((2 * np.asarray(idx, dtype=np.int64) + 1) * ph).astype(np.int64))
    return pc, ph

def modulus_excl(lo, hi, pc, ph):
    lo = np.ascontiguousarray(lo, dtype=np.int64); hi = np.ascontiguousarray(hi, dtype=np.int64)
    out = np.zeros(len(lo), dtype=np.uint8)
    lib().modulus_test(len(lo), _p(lo), _p(hi), _p(pc), ph, _p(out))
    return out.astype(bool)

def modulus_excl_idx(idx, L, pc, ph):
    idx = np.ascontiguousarray(idx, dtype=np.int64).reshape(-1, 5)
    out = np.zeros(len(idx), dtype=np.uint8)
    lib().modulus_test_idx(len(idx), _p(idx), L, _p(pc), ph, _p(out))
    return out.astype(bool)

def krawczyk_excl(lo, hi, pc, ph, dep):
    lo = np.ascontiguousarray(lo, dtype=np.int64); hi = np.ascontiguousarray(hi, dtype=np.int64)
    out = np.zeros(len(lo), dtype=np.uint8)
    d = np.ascontiguousarray(np.array(list(dep) + [0, 0], dtype=np.int32)[:max(len(dep), 1)])
    lib().krawczyk_test(len(lo), _p(lo), _p(hi), _p(pc), ph, len(dep), _p(d), _p(out))
    return out.astype(bool)

def mv_excl(lo, hi, pc, ph):
    lo = np.ascontiguousarray(lo, dtype=np.int64); hi = np.ascontiguousarray(hi, dtype=np.int64)
    out = np.zeros(len(lo), dtype=np.uint8)
    lib().mv_test(len(lo), _p(lo), _p(hi), _p(pc), ph, _p(out))
    return out.astype(bool)

def contained(lo, hi, hlo, hhi):
    lo = np.ascontiguousarray(lo, dtype=np.int64); hi = np.ascontiguousarray(hi, dtype=np.int64)
    hlo = np.ascontiguousarray(hlo, dtype=np.int64).reshape(-1, 5); hhi = np.ascontiguousarray(hhi, dtype=np.int64).reshape(-1, 5)
    out = np.zeros(len(lo), dtype=np.uint8)
    lib().contained(len(lo), _p(lo), _p(hi), len(hlo), _p(hlo), _p(hhi), _p(out))
    return out.astype(bool)

def pair_graphs(hlo, hhi):
    hlo = np.ascontiguousarray(hlo, dtype=np.int64).reshape(-1, 5); hhi = np.ascontiguousarray(hhi, dtype=np.int64).reshape(-1, 5)
    n = len(hlo)
    O = np.zeros((n, n), dtype=np.uint8); B = np.zeros((n, n), dtype=np.uint8); s = np.zeros(n, dtype=np.uint8)
    lib().pair_bounds(n, _p(hlo), _p(hhi), _p(O), _p(B), _p(s))
    return O.astype(bool), B.astype(bool), s.astype(bool)

def choose_dep(fam, pc):
    """dependent parameter pair for the Krawczyk test (numerics only: any choice is valid)"""
    if fam.nt == 0: return ()
    p = pc * (2 * np.pi / NG)
    terms = np.exp(1j * (fam.CE @ p))
    D = np.array([1j * (fam.CE[:, i] * terms).sum() for i in range(fam.m)])
    best, bd = (0, 1), -1
    for a in range(fam.m):
        for b in range(a + 1, fam.m):
            d = abs(D[a].real * D[b].imag - D[a].imag * D[b].real)
            if d > bd: bd, best = d, (a, b)
    return best

# ------------------------------------------------------------------------------------------------ (2) selector
def ordering_ok(i1, i2, lev):
    """exact integer test (necessary condition): some t1 in [i1,i1+1), t2 in [i2,i2+1) (units 2pi/2^lev) with
    t1 <= t2 <= (-t1-t2) mod 2^lev."""
    N = 2 ** lev
    i1 = np.asarray(i1, dtype=np.int64); i2 = np.asarray(i2, dtype=np.int64)
    top = -(i1 + i2)                     # -(t1 + t2) ranges over (top - 2, top]; its integers are top - 1, top
    r = np.mod(top, N)                   # numpy mod: result in [0, N)
    wraps = (r == 0) | (r == 1)          # a multiple of N lies in (top-2, top]: the residue approaches N
    sup = np.where(wraps, N, r)
    return (i1 <= i2) & (sup >= i2)

_S3 = None
def _sqrt3_bounds():
    global _S3
    if _S3 is None:
        import mpmath as mp
        s3lo = float(np.nextafter(np.sqrt(3.0), -np.inf)); s3hi = float(np.nextafter(np.sqrt(3.0), np.inf))
        with mp.workdps(60):
            if not mp.mpf(s3lo) < mp.sqrt(3) < mp.mpf(s3hi):
                raise RuntimeError('sqrt3 enclosure failed')
        _S3 = (s3lo, s3hi)
    return _S3

def select_szollosi(idx, lev):
    idx = np.ascontiguousarray(idx, dtype=np.int64).reshape(-1, 4)
    eC = np.zeros(len(idx), dtype=np.uint8); eS = np.zeros(len(idx), dtype=np.uint8)
    s3lo, s3hi = _sqrt3_bounds()
    lib().szollosi_param_test(len(idx), _p(idx), lev, s3lo, s3hi, _p(eC), _p(eS))
    oa = ordering_ok(idx[:, 0], idx[:, 1], lev); ob = ordering_ok(idx[:, 2], idx[:, 3], lev)
    return oa & ob & (eC == 0) & (eS == 0)

def select_fourier(idx, lev):
    """exact: the box [i,i+1) x [j,j+1) (units 1/2^lev of (a,b)) can meet T = {0<=a<=1/6, 0<=b<=a/2}
    only if i/N <= 1/6 and j/N <= min((i+1)/N, 1/6)/2 ; we keep 6i <= N and 2j <= i+1 (superset)."""
    idx = np.asarray(idx, dtype=np.int64); N = 2 ** lev
    return (6 * idx[:, 0] <= N) & (2 * idx[:, 1] <= idx[:, 0] + 1)

def selector(fam):
    return select_szollosi if fam.name == 'szollosi' else select_fourier

def selected_boxes(fam, level, level_start=3):
    m = fam.m; sel = selector(fam)
    cur = np.indices((2 ** level_start,) * m).reshape(m, -1).T.astype(np.int64)
    cur = cur[sel(cur, level_start)]
    offs = np.indices((2,) * m).reshape(m, -1).T
    for lev in range(level_start + 1, level + 1):
        ch = (cur[:, None, :] * 2 + offs[None]).reshape(-1, m)
        cur = ch[sel(ch, lev)]
    return cur

def prove_empty(fam, idx, lev, max_lev=14):
    """True if the checker's selector excludes every descendant of the box (idx, lev) by level max_lev"""
    sel = selector(fam); m = fam.m
    cur = np.asarray(idx, dtype=np.int64)[None, :]
    offs = np.indices((2,) * m).reshape(m, -1).T
    for lev2 in range(lev + 1, max_lev + 1):
        ch = (cur[:, None, :] * 2 + offs[None]).reshape(-1, m)
        cur = ch[sel(ch, lev2)]
        if len(cur) == 0: return True, lev2
        if len(cur) > 200000: return False, lev2
    return False, max_lev

# ------------------------------------------------------------------------------------------------ (3) coverage
_O5 = np.indices((2,) * 5).reshape(5, -1).T

def children(idx):
    return (idx[:, None, :] * 2 + _O5[None]).reshape(-1, 5)

def dyadic(idx, L):
    s = 2 ** (16 - L)
    return idx * s, (idx + 1) * s

def root_pass(pc, ph, L0=4, L1=7, chunk=2 ** 15):
    """modulus-only branch and bound over the whole phase torus (levels L0..L1) for the parameter box (pc, ph).
    Children are generated and tested chunk by chunk (memory proportional to the survivors).  Returns
    (levels, counts) with levels[L] = (idx_L, parent_L): the surviving boxes of level L and, for L > L0, the row
    of each box's parent in idx_{L-1}."""
    idx = np.indices((2 ** L0,) * 5).reshape(5, -1).T.astype(np.int64)
    idx = idx[~modulus_excl_idx(idx, L0, pc, ph)]
    levels = {L0: (idx, None)}
    counts = [(L0, 2 ** (5 * L0), int(len(idx)))]
    for L in range(L0 + 1, L1 + 1):
        parts = []; pars = []; ntested = 0
        for s in range(0, len(idx), chunk):
            ch = children(idx[s:s + chunk]); ntested += len(ch)
            keep = ~modulus_excl_idx(ch, L, pc, ph)
            parts.append(ch[keep]); pars.append(s + np.nonzero(keep)[0] // 32)
        idx = np.concatenate(parts) if parts else np.zeros((0, 5), dtype=np.int64)
        par = np.concatenate(pars) if pars else np.zeros(0, dtype=np.int64)
        levels[L] = (idx, par)
        counts.append((L, ntested, int(len(idx))))
    return levels, counts

def descend_modulus(levels, pc, ph, Lstart=5, Lend=7):
    """for a parameter sub-box (pc, ph) of the root: the level-Lend root survivors that survive the modulus test
    for (pc, ph) at every level Lstart..Lend (a box excluded at a coarse level has all sub-boxes excluded)"""
    alive = None
    for L in range(Lstart, Lend + 1):
        idx, par = levels[L]
        cand = np.ones(len(idx), dtype=bool) if alive is None else alive[par]
        sub = np.nonzero(cand)[0]
        ok = ~modulus_excl_idx(idx[sub], L, pc, ph)
        alive = np.zeros(len(idx), dtype=bool); alive[sub[ok]] = True
    return levels[Lend][0][alive]

def verify_coverage(fam, idx, lev, hlo, hhi, levels, L7=7, Lcap=12):
    """True iff every phase box that may contain a common unbiased vector for some p in the parameter box
    (idx, lev) is contained in a hull (hlo, hhi).  `levels` is the root pass of an ancestor parameter box.
    Returns (ok, stats)."""
    pc, ph = pbox(idx, lev)
    dep = choose_dep(fam, pc)
    X = descend_modulus(levels, pc, ph, L7, L7); L = L7          # level-L7 root survivors passing the modulus test
    stats = {'dep': list(dep)}
    hlo = np.asarray(hlo, dtype=np.int64).reshape(-1, 5); hhi = np.asarray(hhi, dtype=np.int64).reshape(-1, 5)
    first = True
    while True:
        parts = []; n_in = n_mod = n_unc = 0
        for s in range(0, len(X), 2 ** 14):
            Y = X[s:s + 2 ** 14] if first else children(X[s:s + 2 ** 14])
            n_in += len(Y)
            if not first:                                  # (the first level was modulus-tested by descend_modulus)
                Y = Y[~modulus_excl_idx(Y, L, pc, ph)]
            lo, hi = dyadic(Y, L); n_mod += len(Y)
            if len(hlo) and len(Y):
                c = contained(lo, hi, hlo, hhi)
                Y = Y[~c]; lo = lo[~c]; hi = hi[~c]
            n_unc += len(Y)
            if len(Y):
                v = mv_excl(lo, hi, pc, ph)
                Y = Y[~v]; lo = lo[~v]; hi = hi[~v]
            if len(Y):
                k = krawczyk_excl(lo, hi, pc, ph, dep)
                Y = Y[~k]
            parts.append(Y)
        X = np.concatenate(parts) if parts else np.zeros((0, 5), dtype=np.int64)
        first = False
        stats[f'L{L}'] = (n_in, n_mod, n_unc, int(len(X)))
        if len(X) == 0: return True, stats
        if L >= Lcap: stats['fail_boxes'] = X[:10].tolist(); return False, stats
        if len(X) > 300000: stats['fail'] = 'too many boxes'; return False, stats
        L += 1

# ------------------------------------------------------------------------------------------------ (4) graphs
def cliques6(adj):
    """all 6-cliques (increasing vertex order) of the graph with boolean adjacency adj"""
    n = adj.shape[0]
    nb = [set(np.nonzero(adj[i])[0].tolist()) for i in range(n)]
    out = []
    def ext(cl, cand):
        if len(cl) == 6: out.append(tuple(cl)); return
        for v in sorted(cand):
            if len(cl) + 1 + sum(1 for w in cand if w > v) < 6: break
            ext(cl + [v], {w for w in cand if w > v and w in nb[v]})
    ext([], set(range(n)))
    return out, nb

def verify_graphs(hlo, hhi):
    """True iff every hull passes the self-check and no two 6-cliques of O are joined across B"""
    if len(hlo) == 0: return True, {'n': 0}
    O, B, selfok = pair_graphs(hlo, hhi)
    if not selfok.all(): return False, {'reason': 'self-check', 'n': len(hlo)}
    cl, nbO = cliques6(O)
    nbB = [set(np.nonzero(B[i])[0].tolist()) for i in range(len(hlo))]
    clset = set(cl)
    for K1 in cl:
        common = set.intersection(*[nbB[v] for v in K1])
        if len(common) < 6: continue
        # any 6-clique of O inside `common`?
        cand = sorted(common)
        sub, _ = cliques6(O[np.ix_(cand, cand)])
        if sub: return False, {'reason': 'quartet candidate', 'K1': K1, 'K2': [cand[i] for i in sub[0]], 'n': len(hlo)}
    return True, {'n': len(hlo), 'cliques': len(cl), 'edgesO': int(O.sum() // 2), 'edgesB': int(B.sum() // 2)}
