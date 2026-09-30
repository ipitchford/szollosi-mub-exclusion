"""
selector.py -- rigorous parameter-box selection (GENERATOR side; the checker has its own implementation).
Referee comment 4.2: no library trigonometry, no floating-point angular constants; ordering tests are exact
integer arithmetic; (C) and (S) use table values with explicit error bounds and one-ulp outward steps.

Szöllősi fundamental domain F (paper, Prop. 6): (C) x+y+z+u+v+w = 0, (S) 0 <= arg(x+y+z) <= pi/6 or alpha = 0,
(OA) th_x <= th_y <= th_z, (OB) th_u <= th_v <= th_w, all angles in [0, 2 pi), z = conj(xy), w = conj(uv).
A dyadic box of level l and index i (angles in [i, i+1) 2pi/2^l) is DISCARDED only if one of these conditions
provably fails everywhere on the (closed) box; every point of F lies in the kept box that contains it.
"""
import numpy as np
from rnd import U, TAB_ERR, up, nu, nd
from core import unit_at, NG

S3_DN = float(np.nextafter(np.sqrt(3.0), -np.inf)); S3_UP = float(np.nextafter(np.sqrt(3.0), np.inf))

def ordering_possible_exact(i1, i2, lev):
    """exact integer test: exists t1 in [i1, i1+1), t2 in [i2, i2+1) (units 2pi/2^lev) with t1 <= t2 <= t3,
    t3 = (-t1 - t2) mod N.  Necessary conditions: i1 <= i2 and sup t3 >= i2 (see paper, Sec. 3.4)."""
    N = 2 ** lev
    i1 = np.asarray(i1, dtype=np.int64); i2 = np.asarray(i2, dtype=np.int64)
    b = -(i1 + i2)                      # -(t1+t2) ranges over (b - 2, b]
    # (b-2, b] contains a multiple of N  <=>  floor(b / N) * N > b - 2
    wrap = (np.floor_divide(b, N) * N) > (b - 2)
    sup3 = np.where(wrap, N, b % N)
    return (i1 <= i2) & (sup3 >= i2)

def _dot_lo(n1, n2, v1, v2):
    return nd(nd(n1 * v1) + nd(n2 * v2))

def _dot_hi(n1, n2, v1, v2):
    return nu(nu(n1 * v1) + nu(n2 * v2))

def szollosi_test(idx, lev):
    """keep-mask for level-lev boxes idx (N,4) (rigorous, see module docstring)"""
    idx = np.asarray(idx, dtype=np.int64)
    half = 2 ** (16 - lev)
    c = (2 * idx + 1) * half                              # centres on G17
    hw = float(nu(half * (float(np.nextafter(np.pi, np.inf)) / 2 ** 16)))   # >= pi / 2^lev
    x, y, u, v = [unit_at(c[:, k]) for k in range(4)]
    z = unit_at(-c[:, 0] - c[:, 1]); w = unit_at(-c[:, 2] - c[:, 3])
    terms = [x, y, z, u, v, w]
    S = x + y + z + u + v + w
    Sr = up(6 * TAB_ERR + 6 * U * sum(np.abs(t) for t in terms))
    D = [1j * (x - z), 1j * (y - z), 1j * (u - w), 1j * (v - w)]
    Dr = up(2 * TAB_ERR + 4 * U * 2.0)                      # |x|,|z| <= 1 + TAB_ERR
    # ---- (C): support-function exclusion with candidate directions (any direction is valid)
    gens = np.stack([np.stack([d.real, d.imag], axis=1) for d in D], axis=1) * hw     # numerics only
    P = np.stack([-S.real, -S.imag], axis=1)
    cands = []
    for k in range(4):
        g = gens[:, k]; cands.append(np.stack([-g[:, 1], g[:, 0]], axis=1)); cands.append(np.stack([g[:, 1], -g[:, 0]], axis=1))
    # vertex-to-point directions of the zonotope (numerics only)
    gg = gens.copy(); flip = gg[:, :, 1] < 0; gg[flip] *= -1
    order = np.argsort(np.arctan2(gg[:, :, 1], gg[:, :, 0]), axis=1)
    gg = np.take_along_axis(gg, order[:, :, None], axis=1)
    p = -gg.sum(axis=1)
    for i in range(8):
        cands.append(P - p)
        p = p + 2 * gg[:, i] if i < 4 else p - 2 * gg[:, i - 4]
    rem = nu(6.0 * nu(hw * hw))                               # Taylor remainder 6 hw^2
    exclC = np.zeros(len(idx), dtype=bool)
    for n in cands:
        n1 = n[:, 0]; n2 = n[:, 1]
        nn = nu(np.sqrt(nu(nu(n1 * n1) + nu(n2 * n2))))      # >= |n|
        num = _dot_lo(n1, n2, -S.real, -S.imag)               # <n, -S_fl>
        num = nd(num - nu(nn * Sr))                            # <n, -S_exact> >= ...
        h = np.zeros(len(idx))
        for d in D:
            lo_ = _dot_lo(n1, n2, d.real, d.imag); hi_ = _dot_hi(n1, n2, d.real, d.imag)
            h = nu(h + nu(np.maximum(np.abs(lo_), np.abs(hi_)) + nu(nn * Dr)))
        num = nd(num - nu(hw * h))
        exclC |= (num > nu(nn * rem)) & (nn > 0)
    # ---- (S): alpha over the box lies in the disc (alpha_fl, r), r = 4 hw + error
    al = x + y + z
    r = up(4 * hw + 3 * TAB_ERR + 3 * U * 3.0)
    exclS1 = nu(al.imag + r) < 0                                # Im alpha < 0 everywhere
    s3im = np.where(al.imag >= 0, nd(S3_DN * al.imag), nd(S3_UP * al.imag))
    val = nd(s3im - al.real)                                   # <= sqrt3 Im - Re (float alpha)
    exclS2 = val > nu(2 * r)                                   # (sqrt3 Im - Re)/2 > r
    exclS3 = nu(al.real + r) < 0                                # Re alpha < 0 everywhere (sector lies in Re >= 0)
    # ---- (OA), (OB) exact
    okA = ordering_possible_exact(idx[:, 0], idx[:, 1], lev)
    okB = ordering_possible_exact(idx[:, 2], idx[:, 3], lev)
    return okA & okB & ~exclC & ~exclS1 & ~exclS2 & ~exclS3

def fourier_test(idx, lev):
    """exact rational test: the box [i, i+1) x [j, j+1) (units 1/2^lev of a, b) may meet the triangle
    T = {0 <= a <= 1/6, 0 <= b <= a/2}: kept iff 6 i <= N and 2 j <= i + 1 (a superset of the exact condition)."""
    idx = np.asarray(idx, dtype=np.int64); N = 2 ** lev
    i = idx[:, 0]; j = idx[:, 1]
    return (6 * i <= N) & (2 * j <= i + 1)

def select(test, m, level_max, level_start=3):
    n = 2 ** level_start
    cur = np.indices((n,) * m).reshape(m, -1).T.astype(np.int64)
    cur = cur[test(cur, level_start)]
    offs = np.indices((2,) * m).reshape(m, -1).T
    counts = [(level_start, len(cur))]
    for lev in range(level_start + 1, level_max + 1):
        ch = (cur[:, None, :] * 2 + offs[None, :, :]).reshape(-1, m)
        cur = ch[test(ch, lev)]
        counts.append((lev, len(cur)))
    return cur, counts
