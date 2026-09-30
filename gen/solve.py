"""
solve.py -- per-parameter-box certificate (GENERATOR).

Rigour: a phase box is discarded only by core.modulus_exclude / core.mv_modulus_exclude / core.krawczyk, each
of which proves that the box contains no unbiased vector for any parameter of P (on the constraint surface).
Numerical solutions (Newton) are used ONLY to decide which surviving boxes are grouped into which hull.
Every surviving box of the final level is contained in some hull (hulls.cyclic_hull, repaired), so the hulls
cover U(H_p) for every p in P.  No contraction step is used.  The hulls are returned for storage so that the
independent checker (../chk) can re-verify coverage and the combinatorial step.
"""
import numpy as np, time
import core, hulls, comb
from core import dyadic_boxes, box_geometry, half_widths, modulus_exclude, mv_modulus_exclude, krawczyk, NG16

TWO_PI = 2 * np.pi
KCHUNK = 2 ** 14
_OFFS5 = np.indices((2,) * 5).reshape(5, -1).T

def children_idx(idx):
    return (idx[:, None, :] * 2 + _OFFS5[None, :, :]).reshape(-1, 5)

def cyc_dist(a, b):
    d = np.abs(((a - b) + np.pi) % TWO_PI - np.pi)
    return d.max(axis=-1)

# --------------------------------------------------------------------------------------------- numerics only
def on_surface_point(P, iters=30):
    fam = P.fam
    th = (P.cidx * (TWO_PI / 2 ** 17)).astype(float)
    if fam.CE is None:
        return th % TWO_PI
    d = list(P.dep)
    for _ in range(iters):
        C = fam.C(th)
        if abs(C) < 1e-14: break
        terms = np.exp(1j * (fam.CE @ th))
        D = np.array([1j * (fam.CE[:, i] * terms).sum() for i in range(fam.m)])
        M = np.array([[D[d[0]].real, D[d[1]].real], [D[d[0]].imag, D[d[1]].imag]])
        try: step = np.linalg.solve(M, -np.array([C.real, C.imag]))
        except np.linalg.LinAlgError: break
        th[d[0]] += step[0]; th[d[1]] += step[1]
    return th % TWO_PI

def newton_solutions(H, starts, iters=30, tol=1e-11, dedup=1e-7):
    Hs = H.conj().T
    T = np.array(starts, dtype=float).reshape(-1, 5)
    for it in range(iters):
        V = np.exp(1j * np.concatenate([np.zeros((T.shape[0], 1)), T], axis=1))
        W = V @ Hs.T
        G = (np.abs(W) ** 2 - 6.0)[:, 1:]
        J = 2 * np.real(np.conj(W)[:, 1:, None] * Hs[None, 1:, 1:] * (1j * V[:, None, 1:]))
        with np.errstate(all='ignore'):
            try:
                dT = np.linalg.solve(J + 1e-14 * np.eye(5)[None], G[:, :, None])[:, :, 0]
            except np.linalg.LinAlgError:
                dT = np.zeros_like(T)
        dT = np.nan_to_num(dT)
        step = np.clip(np.abs(dT).max(axis=1, keepdims=True), 1e-300, None)
        damp = np.minimum(1.0, 0.5 / step)
        T = (T - damp * dT) % TWO_PI
    V = np.exp(1j * np.concatenate([np.zeros((T.shape[0], 1)), T], axis=1))
    res = np.abs(np.abs(V @ Hs.T) ** 2 - 6.0).max(axis=1)
    T = T[res < tol]
    out = []
    for t in T:
        if not out or np.abs(np.exp(1j * t) - np.exp(1j * np.array(out))).max(axis=1).min() > dedup:
            out.append(t)
    return np.array(out).reshape(-1, 5)

def multistart_solutions(H, nstart=20000, seed=0):
    rng = np.random.default_rng(seed)
    return newton_solutions(H, rng.uniform(0, TWO_PI, (nstart, 5)))

# --------------------------------------------------------------------------------------------- rigorous B&B
def test_level(idx, L, P, use_mv=True, use_k=True, chunk=2 ** 18):
    """keep-mask of dyadic phase boxes (level L) not excluded by the rigorous tests"""
    keep = np.zeros(len(idx), dtype=bool)
    for s in range(0, len(idx), chunk):
        lo, hi = dyadic_boxes(idx[s:s + chunk], L)
        V, rho = box_geometry(lo, hi)
        k = ~modulus_exclude(V, rho, P)
        sub = np.nonzero(k)[0]
        if use_mv and len(sub):
            k2 = ~mv_modulus_exclude(V[sub], rho[sub], P); k[sub[~k2]] = False
        sub = np.nonzero(k)[0]
        if use_k and len(sub):
            hw5 = half_widths(lo[sub], hi[sub])
            for s2 in range(0, len(sub), KCHUNK):
                ss = sub[s2:s2 + KCHUNK]; hh = np.arange(s2, min(s2 + KCHUNK, len(sub)))
                ex = krawczyk(V[ss], rho[ss], hw5[hh], P)
                k[ss[ex]] = False
        keep[s:s + chunk] = k
    return keep

def root_pass(P, L0=4, L1=7, chunk=2 ** 19):
    """modulus-only B&B over the full phase torus from level L0 to L1 for the (root) parameter box P"""
    n0 = 2 ** L0
    idx = np.indices((n0,) * 5).reshape(5, -1).T.astype(np.int64)
    L = L0
    counts = []
    while True:
        keep = test_level(idx, L, P, use_mv=False, use_k=False, chunk=chunk)
        counts.append((L, len(idx), int(keep.sum())))
        idx = idx[keep]
        if L >= L1: break
        idx = children_idx(idx); L += 1
    return idx, L, counts

def _declared_radii(sols, P, Lmax, Rfac=1.3, Rmin=0.02, Rcap=0.15):
    if len(sols) == 0: return np.zeros(0)
    nd_ = 2 ** Lmax
    sidx = np.floor(sols / (TWO_PI / nd_)).astype(np.int64) % nd_
    lo, hi = dyadic_boxes(sidx, Lmax); V, rho = box_geometry(lo, hi)
    ex, Ks, Kr, ok = krawczyk(V, rho, half_widths(lo, hi), P, return_K=True)
    R = np.maximum(Rfac * (Kr[:, :5].max(axis=1) + np.abs(Ks[:, :5]).max(axis=1)) + 2 * np.pi / nd_, Rmin)
    R = np.where(ok & np.isfinite(R), R, 0.3)
    return np.minimum(R, Rcap)

def _cover(idx, L, sols, R, covered):
    if len(idx) == 0 or len(sols) == 0: return idx
    c = (idx + 0.5) * (TWO_PI / 2 ** L)
    d = cyc_dist(c[:, None, :], sols[None, :, :])
    within = d <= R[None, :]
    cov = within.any(axis=1); owner = np.argmax(within, axis=1)
    for i in range(len(sols)):
        m = cov & (owner == i)
        if m.any():
            lo, hi = dyadic_boxes(idx[m], L); covered[i].append((lo, hi))
    return idx[~cov]

def child_box(P, idx7, L7, sols_seed, Lmax=9, Ld=7, max_boxes=2_000_000, n_rescue=3000, rescue_rounds=3):
    t0 = time.time()
    info = {'lev': P.lev, 'idx': P.idx.tolist()}
    keep = test_level(idx7, L7, P)
    idx = idx7[keep]
    info['L7'] = (int(len(idx7)), int(len(idx)))
    th = on_surface_point(P)
    H = P.fam.H(th)
    sols = newton_solutions(H, sols_seed) if len(sols_seed) else np.zeros((0, 5))
    if len(sols) < 24:
        sols = multistart_solutions(H, nstart=20000, seed=1)
    R = _declared_radii(sols, P, Lmax)
    L = L7
    covered = [[] for _ in range(len(sols))]
    while True:
        if L > L7:
            keep = test_level(idx, L, P); idx = idx[keep]
        if L >= Ld:
            idx = _cover(idx, L, sols, R, covered)
        if L >= Lmax or len(idx) == 0: break
        if len(idx) * 32 > max_boxes:
            info.update({'ok': False, 'reason': 'too many boxes', 'time': round(time.time() - t0, 1)}); return info
        idx = children_idx(idx); L += 1
    n_rescued = 0
    for _ in range(rescue_rounds):
        if len(idx) == 0: break
        c = (idx + 0.5) * (TWO_PI / 2 ** L)
        if len(c) > n_rescue:
            c = c[np.random.default_rng(len(idx)).choice(len(c), n_rescue, replace=False)]
        new = newton_solutions(H, c, iters=60)
        if len(new) and len(sols):
            new = new[cyc_dist(new[:, None, :], sols[None, :, :]).min(axis=1) > 1e-6]
        if len(new) == 0: break
        Rn = _declared_radii(new, P, Lmax)
        sols = np.concatenate([sols, new]) if len(sols) else new
        R = np.concatenate([R, Rn]); covered.extend([[] for _ in range(len(new))])
        n_rescued += len(new)
        idx = _cover(idx, L, sols, R, covered)
    info['nsol'] = int(len(sols)); info['rescued'] = n_rescued; info['leftover'] = int(len(idx))
    hl_list = []; hh_list = []
    for i in range(len(sols)):
        if not covered[i]: continue                     # a hull must contain surviving boxes to be needed
        lo = np.concatenate([a for a, b in covered[i]]); hi = np.concatenate([b for a, b in covered[i]])
        hl, hh = hulls.cyclic_hull(lo, hi); hl_list.append(hl); hh_list.append(hh)
    if len(idx):
        lo, hi = dyadic_boxes(idx, L)
        Lc = max(L - 2, 4)
        cidx = ((lo + hi) // 2 // (2 ** (16 - Lc))) % (2 ** Lc)
        comps = hulls.cluster(cidx, Lc)
        for c in comps:
            hl, hh = hulls.cyclic_hull(lo[c], hi[c]); hl_list.append(hl); hh_list.append(hh)
        info['extra_hulls'] = len(comps)
    hl = np.array(hl_list, dtype=np.int64).reshape(-1, 5); hh = np.array(hh_list, dtype=np.int64).reshape(-1, 5)
    info['hulls'] = int(len(hl))
    info['hull_lo'] = hl.tolist(); info['hull_hi'] = hh.tolist()
    if len(hl) == 0:
        info.update({'ok': True, 'cert': {'n': 0, 'cliques': 0}, 'time': round(time.time() - t0, 1)})
        info['sols'] = sols; return info
    V, rho = box_geometry(hl, hh)
    r5 = half_widths(hl, hh)
    info['max_hw'] = float(r5.max())
    if not comb.self_ok(r5):
        info.update({'ok': False, 'reason': 'hull too wide', 'time': round(time.time() - t0, 1)}); info['sols'] = sols; return info
    ok, cinfo = comb.certificate(V, rho)
    info.update({'ok': ok, 'cert': cinfo, 'time': round(time.time() - t0, 1)})
    info['sols'] = sols
    return info
