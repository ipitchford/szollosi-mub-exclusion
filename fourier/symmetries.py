"""
symmetries.py -- exact search for parameter symmetries of the Fourier family F(a,b) (Jaming et al. 2009, eq. 3).

Every entry of sqrt6 F(a,b) is eta^z x^i y^j with eta = exp(i pi/6), x = exp(2 pi i a), y = exp(2 pi i b); we
represent it by the exponent triple (z mod 12, i, j).  An affine map T(a,b) = (p a + q b + m/12, r a + t b + n/12)
with integer p,q,r,t is a SYMMETRY if for all (a,b)
        F(T(a,b)) = D1 P1 G(F(a,b)) P2 D2                                               (*)
where P1, P2 are permutation matrices, D1, D2 diagonal unitary (functions of (a,b)) and G is one of
id, entrywise conjugation, transposition, conjugate transposition.  Each such G, and left/right multiplication by
monomial unitaries, maps a MUB quartet containing (Id, F) to one containing (Id, F'), so symmetries preserve
extendability to a quartet.  In exponent notation (*) says that
        R_kj = E'_kj - E_{pi(k), sigma(j)}   (in Z12 x Z x Z)
is of the form alpha_k + beta_j; the check below is exact integer arithmetic.
"""
import numpy as np, itertools, json, sys
from fractions import Fraction as Fr

ROWS = [["1", "1", "1", "1", "1", "1"],
        ["1", "-w2 x", "w", "-x", "w2", "-w x"],
        ["1", "w y", "w2", "y", "w", "w2 y"],
        ["1", "-1", "1", "-1", "1", "-1"],
        ["1", "w2 x", "w", "x", "w2", "w x"],
        ["1", "-w y", "w2", "-y", "w", "-w2 y"]]

def exponents():
    E = np.zeros((6, 6, 3), dtype=np.int64)       # (z mod 12, i, j): eta^z x^i y^j, omega = eta^4, -1 = eta^6
    for k in range(6):
        for j in range(6):
            s = ROWS[k][j]; z = 0; i = jj = 0
            if s.startswith('-'): z += 6; s = s[1:]
            for tk in s.split():
                if tk == 'w': z += 4
                elif tk == 'w2': z += 8
                elif tk == 'x': i += 1
                elif tk == 'y': jj += 1
            E[k, j] = (z % 12, i, jj)
    return E

E0 = exponents()
PERMS = np.array(list(itertools.permutations(range(6))), dtype=np.int64)     # 720 x 6

def transformed(E, T):
    """exponents of F(T(a,b)) in terms of (a,b): T = (p,q,r,t,m,n)"""
    p, q, r, t, m, n = T
    z, i, j = E[..., 0], E[..., 1], E[..., 2]
    return np.stack([(z + i * m + j * n) % 12, i * p + j * r, i * q + j * t], axis=-1)

def G(E, g):
    if g in (1, 3): E = np.stack([(-E[..., 0]) % 12, -E[..., 1], -E[..., 2]], axis=-1)
    if g in (2, 3): E = np.transpose(E, (1, 0, 2))
    return E

def norm(v):
    """normalise exponent vectors: z mod 12"""
    v = v.copy(); v[..., 0] %= 12; return v

def find_equivalence(Et, Es):
    """find (pi, sigma) with Et_kj - Es_{pi(k),sigma(j)} = alpha_k + beta_j; returns (pi, sigma, alpha, beta) or None"""
    # for each pi: A[pi, k, l, :] = Et[k, :, :] vs Es[pi(k), l]
    for pi in PERMS:
        Esp = Es[pi]                                              # rows permuted: Esp[k, l] = Es[pi(k), l]
        # candidate sigma(0) = l0: alpha_k = Et[k,0] - Esp[k,l0] (beta_0 = 0)
        for l0 in range(6):
            alpha = Et[:, 0] - Esp[:, l0]
            sigma = [-1] * 6; used = set(); ok = True
            for j in range(6):
                found = -1
                for l in range(6):
                    if l in used: continue
                    d = Et[:, j] - Esp[:, l] - alpha                 # must be constant in k
                    d[:, 0] %= 12
                    if np.all(d == d[0]): found = l; break
                if found < 0: ok = False; break
                sigma[j] = found; used.add(found)
            if ok:
                beta = [(Et[0, j] - Esp[0, sigma[j]] - alpha[0]) for j in range(6)]
                return pi.tolist(), sigma, alpha, beta
    return None

def fast_candidates(Et, Es):
    """cheap necessary condition: multisets of row exponent-differences must match up to permutation (skip)"""
    return True

def search(lin_range=(-2, -1, 0, 1, 2), trans=range(12)):
    found = []
    mats = [M for M in itertools.product(lin_range, repeat=4) if abs(M[0] * M[3] - M[1] * M[2]) == 1]
    for (p, q, r, t) in mats:
        for m in trans:
            for n in trans:
                T = (p, q, r, t, m, n)
                Et = transformed(E0, T)
                for g in range(4):
                    Es = G(E0, g)
                    # quick necessary test: the multiset of (i,j) parts of rows must match up to permutation and sign-free
                    res = find_equivalence(Et, Es)
                    if res is not None:
                        found.append((T, g, res)); break
    return found

if __name__ == '__main__':
    pass

def haagerup_multiset(E):
    """sorted multiset of the invariants E_ki + E_lj - E_kj - E_li (exponent form; invariant under row/column
    permutations and diagonal scalings, and under transposition)"""
    L = E[:, None, :, None, :] + E[None, :, None, :, :] - E[:, None, None, :, :] - E[None, :, :, None, :]
    z = L[..., 0] % 12; i = L[..., 1]; j = L[..., 2]
    code = (z * 1000 + (i + 400)) * 1000 + (j + 400)
    return np.sort(code.ravel())

def search_fast(lin_range=(-2, -1, 0, 1, 2), trans=range(12), log=print):
    H0 = haagerup_multiset(E0)
    H1 = haagerup_multiset(G(E0, 1))
    found = []; ncand = 0
    mats = [M for M in itertools.product(lin_range, repeat=4) if abs(M[0] * M[3] - M[1] * M[2]) == 1]
    for (p, q, r, t) in mats:
        for m in trans:
            for n in trans:
                T = (p, q, r, t, m, n)
                Et = transformed(E0, T)
                Ht = haagerup_multiset(Et)
                for g in range(4):
                    if not np.array_equal(Ht, H1 if g in (1, 3) else H0): continue
                    ncand += 1
                    res = find_equivalence(Et, G(E0, g))
                    if res is not None:
                        found.append((T, g, res)); break
    log(f'{len(mats)} linear parts x {len(trans)**2} translations: {ncand} passed the invariant filter, {len(found)} symmetries')
    return found
