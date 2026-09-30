"""
comb.py -- combinatorial certificate on hulls (generator side), with rigorous rounding.

Hull i: centre unit numbers c_ik (from the table, |c - exact| <= TAB_ERR), radii rho_ik >= |exp(i th) - c_ik|
(chord bound), angular half-widths r_ik.  For v in hull i, w in hull j:
   (1/6)|<v,w>|  in  [ |G_ij| - R_ij , |G_ij| + R_ij ],  G_ij = (1/6) sum_k conj(c_ik) c_jk,
   R_ij = (1/6) sum_k ((1+TAB_ERR)(rho_ik + rho_jk) + rho_ik rho_jk) + rounding.
For v, w in the same hull i: (1/6)|<v,w>| >= 1 - (1/3) sum_k r_ik^2   (cos t >= 1 - t^2/2).
O (possibly orthogonal): NOT certainly |G| - R > 0.   B (possibly unbiased): NOT certainly ||G| - 1/sqrt6| > R.
"""
import numpy as np
from rnd import U, TAB_ERR, up, nu, nd, INV_SQ6_UP, INV_SQ6_DN

def pair_intervals(C, R):
    n = C.shape[0]
    G = (np.conj(C) @ C.T) / 6.0
    absC = np.abs(C)
    rad = up(((1 + TAB_ERR) * (R @ np.ones((6, n)) + np.ones((n, 6)) @ R.T) + R @ R.T) / 6.0
             + 12 * U * (absC @ absC.T) / 6.0)
    a = np.abs(G)
    rad = up(rad + 4 * U * a)
    return a, rad

def self_ok(r5):
    """r5: (n,5) angular half-widths (upper bounds).  True iff every hull certainly has 1 - sum r^2/3 > 1/sqrt6."""
    s = up((r5 * r5).sum(axis=1) / 3.0)
    return bool((nd(1.0 - s) > INV_SQ6_UP).all())

def graphs(C, R):
    a, rad = pair_intervals(C, R)
    notO = nd(a * (1 - 2 * U)) > rad                         # certainly |G| > 0
    # certainly | |G| - 1/sqrt6 | > rad :  |G| < 1/sqrt6 - rad   or   |G| > 1/sqrt6 + rad
    below = nu(a * (1 + 2 * U)) < nd(INV_SQ6_DN - rad)
    above = nd(a * (1 - 2 * U)) > nu(INV_SQ6_UP + rad)
    notB = below | above
    O = ~notO; B = ~notB
    np.fill_diagonal(O, False); np.fill_diagonal(B, False)
    return O, B

def bits(adj):
    """adjacency rows as Python integers; indices converted to Python int BEFORE shifting (see paper 8.3)"""
    return [sum(1 << int(j) for j in np.nonzero(adj[i])[0]) for i in range(adj.shape[0])]

def cliques6_iter(bitsO, n, restrict=None):
    full = (1 << n) - 1 if restrict is None else restrict
    stack = [([], full)]
    while stack:
        cl, cand = stack.pop()
        if len(cl) == 6:
            yield cl; continue
        c = cand
        while c:
            v = (c & -c).bit_length() - 1
            c &= c - 1
            stack.append((cl + [v], c & bitsO[v]))

def certificate(C, R, max_cliques=5_000_000):
    """True iff no two 6-cliques of O are joined across B (all 36 pairs)."""
    n = C.shape[0]
    if n == 0: return True, {'n': 0, 'cliques': 0}
    O, B = graphs(C, R)
    bO = bits(O); bB = bits(B)
    ncl = 0
    for c in cliques6_iter(bO, n):
        ncl += 1
        if ncl > max_cliques:
            return False, {'reason': 'too many cliques', 'n': n, 'cliques': ncl}
        common = (1 << n) - 1
        for v in c: common &= bB[v]
        if bin(common).count('1') < 6: continue
        for d in cliques6_iter(bO, n, restrict=common):
            return False, {'reason': 'quartet candidate', 'C': c, 'D': d, 'n': n, 'cliques': ncl}
    return True, {'n': n, 'cliques': ncl, 'edgesO': int(O.sum() // 2), 'edgesB': int(B.sum() // 2)}
