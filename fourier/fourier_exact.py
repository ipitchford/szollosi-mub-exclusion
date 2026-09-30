"""fourier_exact.py -- exact (SymPy) check that sqrt6 F(a,b) of Jaming et al. (2009, eq. 3) satisfies F F^* = 6 I
identically in the unimodular symbols x = e^{2 pi i a}, y = e^{2 pi i b}."""
import sympy as sp, json
x, y = sp.symbols('x y', nonzero=True)
w = sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2
F = sp.Matrix([[1, 1, 1, 1, 1, 1],
               [1, -w**2 * x, w, -x, w**2, -w * x],
               [1, w * y, w**2, y, w, w**2 * y],
               [1, -1, 1, -1, 1, -1],
               [1, w**2 * x, w, x, w**2, w * x],
               [1, -w * y, w**2, -y, w, -w**2 * y]])
conj = lambda e: sp.conjugate(e).subs({sp.conjugate(x): 1 / x, sp.conjugate(y): 1 / y})
G = sp.simplify(F * F.T.applyfunc(conj) - 6 * sp.eye(6))
ok = all(sp.simplify(sp.expand(e)) == 0 for e in G)
unimod = all(sp.simplify(sp.expand(e * conj(e))) == 1 for e in F)
print(json.dumps({'F_Fstar_equals_6I_identically': ok, 'entries_unimodular': unimod}))
json.dump({'F_Fstar_equals_6I_identically': ok, 'entries_unimodular': unimod}, open('fourier_exact.json', 'w'))
