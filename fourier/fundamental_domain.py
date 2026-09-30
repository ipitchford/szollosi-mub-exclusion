"""
fundamental_domain.py -- the triangle T = conv{(0,0), (1/6,0), (1/6,1/12)} is a fundamental domain for the
parameter symmetries of the Fourier family (exact rational arithmetic).

1. search: all affine maps (a,b) -> (p a + q b + m/12, r a + t b + n/12), |p|,|q|,|r|,|t| <= 2, det = +-1, that
   are symmetries of F in the sense of symmetries.py (exact exponent arithmetic).
2. independent numerical confirmation of every symmetry at random parameters (numpy, residual ~1e-15).
3. group structure: closure under composition and inverses (mod Z^2).
4. covering: the preimages T_g^{-1}(T) have total area 1 and pairwise disjoint interiors on the torus
   (separating-axis test in exact rationals), hence cover the torus (a closed set of full measure).
   Consequently every (a,b) has a symmetry image in T, and a MUB quartet containing (Id, F(a,b)) exists for
   some (a,b) iff one exists for some (a,b) in T.
"""
import numpy as np, itertools, json, time, sys
from fractions import Fraction as Fr
from symmetries import E0, search_fast, transformed, G

def numeric_check(T, g, res, trials=20, rng=np.random.default_rng(7)):
    pi, sigma, alpha, beta = res
    eta = np.exp(1j * np.pi / 6)
    def M(E, a, b):
        return eta ** E[..., 0] * np.exp(2j * np.pi * (E[..., 1] * a + E[..., 2] * b))
    worst = 0.0
    for _ in range(trials):
        a, b = rng.uniform(0, 1, 2)
        p, q, r, t, m, n = T
        ap, bp = p * a + q * b + m / 12, r * a + t * b + n / 12
        Fp = M(E0, ap, bp)
        Fs = M(G(E0, g), a, b)
        d1 = np.array([eta ** (al[0] % 12) * np.exp(2j * np.pi * (al[1] * a + al[2] * b)) for al in alpha])
        d2 = np.array([eta ** (be[0] % 12) * np.exp(2j * np.pi * (be[1] * a + be[2] * b)) for be in beta])
        rhs = d1[:, None] * Fs[np.ix_(pi, sigma)] * d2[None, :]
        worst = max(worst, float(np.abs(Fp - rhs).max()))
    return worst

def affine(T):
    p, q, r, t, m, n = T
    return (Fr(p), Fr(q), Fr(r), Fr(t), Fr(m, 12), Fr(n, 12))

def apply(A, pt):
    p, q, r, t, s, u = A
    return (p * pt[0] + q * pt[1] + s, r * pt[0] + t * pt[1] + u)

def inverse(A):
    p, q, r, t, s, u = A
    det = p * t - q * r
    ip, iq, ir, it = t / det, -q / det, -r / det, p / det
    return (ip, iq, ir, it, -(ip * s + iq * u), -(ir * s + it * u))

def compose(A, B):          # A o B
    p, q, r, t, s, u = A; P, Q, R, S, s2, u2 = B
    return (p * P + q * R, p * Q + q * S, r * P + t * R, r * Q + t * S, p * s2 + q * u2 + s, r * s2 + t * u2 + u)

def key(A):
    return (A[0], A[1], A[2], A[3], A[4] % 1, A[5] % 1)

TRI = [(Fr(0), Fr(0)), (Fr(1, 6), Fr(0)), (Fr(1, 6), Fr(1, 12))]

def area(P):
    (x1, y1), (x2, y2), (x3, y3) = P
    return abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2

def separated(P, Q):
    """True iff the closed triangles P, Q admit a weakly separating line parallel to an edge (disjoint interiors)"""
    for poly in (P, Q):
        for i in range(3):
            (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % 3]
            nx, ny = y2 - y1, x1 - x2
            pa = [nx * x + ny * y for x, y in P]; pb = [nx * x + ny * y for x, y in Q]
            if max(pa) <= min(pb) or max(pb) <= min(pa): return True
    return False

if __name__ == '__main__':
    log = {}
    t0 = time.time()
    found = search_fast(log=lambda s: print(s, flush=True))
    log['search'] = {'maps_found': len(found), 'time_s': round(time.time() - t0, 1),
                     'variants_used': sorted(set(g for _, g, _ in found))}
    worst = max(numeric_check(T, g, res) for T, g, res in found)
    log['numeric_confirmation_max_abs_residual'] = worst
    print('numerical confirmation, max |F(T(a,b)) - D1 P1 G(F(a,b)) P2 D2| =', worst)
    maps = [affine(T) for T, g, res in found]
    keys = {key(A) for A in maps}
    closed = all(key(compose(A, B)) in keys for A in maps for B in maps)
    inv_closed = all(key(inverse(A)) in keys for A in maps)
    log['group'] = {'order': len(keys), 'closed_under_composition': closed, 'closed_under_inverse': inv_closed}
    print('group: order', len(keys), 'closed', closed, 'inverses', inv_closed)
    pre = [[apply(inverse(A), v) for v in TRI] for A in maps]
    # reduce modulo Z^2: first vertex in [0,1)^2 (then all coordinates lie in (-1, 2) and shifts |u|,|v| <= 3 suffice)
    import math
    pre = [[(x - math.floor(P[0][0]), y - math.floor(P[0][1])) for x, y in P] for P in pre]
    assert all(-1 < x < 2 and -1 < y < 2 for P in pre for x, y in P)
    total = sum(area(P) for P in pre)
    log['total_area_of_preimages'] = str(total)
    bad = 0; npairs = 0
    shifts = [(Fr(u), Fr(v)) for u in range(-3, 4) for v in range(-3, 4)]
    for i in range(len(pre)):
        for j in range(i, len(pre)):
            for (u, v) in shifts:
                if i == j and u == 0 and v == 0: continue
                Q = [(x + u, y + v) for x, y in pre[j]]
                npairs += 1
                if not separated(pre[i], Q): bad += 1
    log['pairs_checked'] = npairs; log['overlapping_pairs'] = bad
    ok = (total == 1) and bad == 0
    log['fundamental_domain_proved'] = ok
    log['symmetries'] = [{'T': list(T), 'variant': g, 'row_perm': res[0], 'col_perm': res[1],
                          'D1_exponents': [list(map(int, a)) for a in res[2]], 'D2_exponents': [list(map(int, b)) for b in res[3]]}
                         for T, g, res in found]
    print('total area', total, '; overlapping pairs', bad, 'of', npairs, '; fundamental domain proved:', ok)
    json.dump(log, open('fundamental_domain.json', 'w'), indent=1)
