"""
lemma5_exact.py -- exact symbolic verification (SymPy) of the matrix identities of Lemma 5 and of the facts about
the symmetry maps used in Proposition 6.  The unimodular variables are treated as Laurent-polynomial symbols
(conj(x) = 1/x), omega = exp(2 pi i/3) as an exact algebraic number; every identity is checked entrywise by
exact simplification to zero.
"""
import sympy as sp, json, itertools
x, y, u, v = sp.symbols('x y u v', nonzero=True)
w3 = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2          # omega = exp(2 pi i / 3)
assert sp.simplify(w3 ** 3 - 1) == 0

def conj(e):
    """conjugation of a Laurent expression in unimodular symbols with exact complex constants"""
    return sp.conjugate(e).subs({sp.conjugate(x): 1 / x, sp.conjugate(y): 1 / y, sp.conjugate(u): 1 / u, sp.conjugate(v): 1 / v})

def circ(a0, a1, a2):
    return sp.Matrix([[a0, a1, a2], [a2, a0, a1], [a1, a2, a0]])

def H(X, Y, U, V):
    Z = 1 / (X * Y); W = 1 / (U * V)
    A = circ(1, 1 / X, Z); B = circ(1, 1 / U, W)
    Bs = B.T.applyfunc(conj); As = A.T.applyfunc(conj)
    return sp.Matrix(sp.BlockMatrix([[A, B], [Bs, -As]]))

def zero(M):
    return all(sp.simplify(sp.expand(e)) == 0 for e in M)

def blk(M1, M2):
    return sp.Matrix(sp.BlockMatrix([[M1, sp.zeros(3)], [sp.zeros(3), M2]]))

I3 = sp.eye(3)
P23 = sp.Matrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
Pi = circ(0, 1, 0)                                        # cyclic shift
Dw = sp.diag(1, w3, w3 ** 2)
Hp = H(x, y, u, v)
res = {}
# maps in terms of (x, y, u, v); z = 1/(xy), w = 1/(uv)
z = 1 / (x * y); w = 1 / (u * v)
maps = {
    'sigma1': (y, z, u, v),                 # (x,y,z) -> (y,z,x)
    'sigma2': (x, y, v, w),                 # (u,v,w) -> (v,w,u)
    'rho': (1 / z, 1 / y, 1 / w, 1 / v),     # (x,y,z,u,v,w) -> (zb, yb, xb, wb, vb, ub)
    'tau': (1 / z, 1 / y, 1 / u, 1 / v),     # -> (zb, yb, xb, ub, vb, wb)
    'kappa': (1 / x, 1 / y, 1 / u, 1 / v),
    'omega': (w3 * x, w3 * y, w3 * u, w3 * v),
    'varsigma': (u, v, x, y),
}
# 1. tau: H(tau p) = H(p)^T
res['1_tau_transpose'] = zero(H(*maps['tau']) - Hp.T)
# 2. kappa: H(kappa p) = conj(H(p))
res['2_kappa_conjugate'] = zero(H(*maps['kappa']) - Hp.applyfunc(conj))
# 3. rho: H(rho p) = (P23 + P23) H (P23 + P23)
Q = blk(P23, P23)
res['3_rho'] = zero(H(*maps['rho']) - Q * Hp * Q)
# 4. omega: H(omega p) = (Dw + Dw) H (Dw + Dw)^-1
Dd = blk(Dw, Dw)
res['4_omega'] = zero(H(*maps['omega']) - Dd * Hp * Dd.inv())
# 5. varsigma: H(varsigma p) = (I + (-I)) H [[0, I], [I, 0]]
J = sp.Matrix(sp.BlockMatrix([[sp.zeros(3), I3], [I3, sp.zeros(3)]]))
res['5_varsigma'] = zero(H(*maps['varsigma']) - blk(I3, -I3) * Hp * J)
# 6. sigma1, sigma2 with the exact exponent signs determined here
found6 = {}
for s in (1, -1):
    Pa = Pi if s == 1 else Pi.T
    Pb = Pi.T if s == 1 else Pi
    if zero(H(*maps['sigma1']) - blk(x * I3, Pa) * Hp * blk(Pb, (1 / x) * I3)): found6['sigma1'] = s
    if zero(H(*maps['sigma2']) - blk(u * I3, Pa) * Hp * blk((1 / u) * I3, Pb)): found6['sigma2'] = s
res['6_sigma1_sign'] = found6.get('sigma1'); res['6_sigma2_sign'] = found6.get('sigma2')
if 'sigma2' not in found6:
    # search all combinations of shifts for the sigma2 identity
    for a, b in itertools.product([I3, Pi, Pi.T], repeat=2):
        for L, R in [(blk(u * I3, a), blk((1 / u) * I3, b)), (blk(a, u * I3), blk(b, (1 / u) * I3)), (blk(u * I3, a), blk(b, (1 / u) * I3)), (blk(a, u * I3), blk((1 / u) * I3, b))]:
            if zero(H(*maps['sigma2']) - L * Hp * R): res['6_sigma2_form'] = [str(L), str(R)]; break
# the constraint: C(p') in terms of C(p)
C = lambda X, Y, U, V: X + Y + 1 / (X * Y) + U + V + 1 / (U * V)
Cp = C(x, y, u, v)
alpha = x + y + z
for name, m in maps.items():
    Cm = C(*m)
    rel = None
    for lab, target in [('C', Cp), ('conj C', conj(Cp)), ('omega C', w3 * Cp)]:
        if sp.simplify(sp.expand(Cm - target)) == 0: rel = lab
    res[f'constraint_{name}'] = rel
# alpha transformations (alpha' = x' + y' + z'):  on Sigma (u+v+w = -alpha)
alpha_rel = {}
for name, m in maps.items():
    X, Y, U, V = m; ap = X + Y + 1 / (X * Y)
    for lab, target in [('alpha', alpha), ('conj alpha', conj(alpha)), ('omega alpha', w3 * alpha), ('u+v+w (= -alpha on Sigma)', u + v + w)]:
        if sp.simplify(sp.expand(ap - target)) == 0: alpha_rel[name] = lab; break
res['alpha_transformations'] = alpha_rel
# label actions used in Proposition 6
def compose(m2, m1):
    X, Y, U, V = m1
    return tuple(sp.simplify(e.subs({x: X, y: Y, u: U, v: V}, simultaneous=True)) for e in m2)
kt = compose(maps['kappa'], maps['tau'])      # kappa after tau
kr = compose(maps['kappa'], maps['rho'])
res['kappa_tau_reverses_first_triple'] = all(sp.simplify(a - b) == 0 for a, b in zip(kt, (z, y, u, v)))    # (z, y, x, u, v, w)
res['kappa_rho_reverses_both'] = all(sp.simplify(a - b) == 0 for a, b in zip(kr, (z, y, w, v)))            # (z, y, x, w, v, u)
# Hadamard property: H H* = 6 I on Sigma  <=>  off-diagonal blocks vanish; check H H* - 6I = circulant(C-terms)
HH = sp.simplify(Hp * Hp.T.applyfunc(conj))
res['HHstar_diagonal_6'] = all(sp.simplify(HH[i, i] - 6) == 0 for i in range(6))
# every off-diagonal entry of H H* is 0, C(p), or conj(C(p))
offd = set()
for i in range(6):
    for j in range(6):
        if i == j: continue
        e = sp.simplify(HH[i, j])
        offd.add('0' if e == 0 else ('C' if sp.simplify(e - Cp) == 0 else ('conjC' if sp.simplify(e - conj(Cp)) == 0 else 'other')))
res['HHstar_offdiagonal_entries'] = sorted(offd)
print(json.dumps(res, indent=1, default=str))
json.dump(res, open('lemma5_exact.json', 'w'), indent=1, default=str)
