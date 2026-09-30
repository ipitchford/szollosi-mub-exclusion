"""
family.py -- 'admissible families' of complex Hadamard matrices of order 6 (generator side).

A family is given by
    H(p)_{kj} = c_{kj} exp(i <E_{kj}, p>),   p in T^m (angles),   c_{kj} unimodular constants,
optionally subject to a trigonometric constraint  C(p) = sum_t exp(i <e_t, p>) = 0  (two real equations).
The certificate machinery only needs (m, E, c, and the constraint exponents); see paper, Section 4.

Two families are defined:
  szollosi : H(p) = [A B; B* -A*], A = circ(1, conj x, conj(xy)), B = circ(1, conj u, conj(uv)),
             p = (th_x, th_y, th_u, th_v), constraint x + y + conj(xy) + u + v + conj(uv) = 0.
  fourier  : F(a,b) of Jaming et al. (2009), eq. (3), times sqrt 6, p = (2 pi a, 2 pi b), no constraint.
The coefficients of the Fourier family are sixth roots of unity; their floating-point values carry an
explicit error bound COEF_ERR (verified with mpmath interval arithmetic at import).
"""
import numpy as np

class Family:
    def __init__(self, name, m, E, coef, coef_err, CE=None):
        self.name = name; self.m = m
        self.E = np.asarray(E, dtype=np.int64)              # (6,6,m)
        self.coef = np.asarray(coef, dtype=complex)          # (6,6)
        self.coef_err = np.asarray(coef_err, dtype=float)    # (6,6) error bounds |coef_fl - coef| <= err
        self.CE = None if CE is None else np.asarray(CE, dtype=np.int64)   # (T, m)
        self.nc = 0 if CE is None else 2
        self.E1 = np.abs(self.E).sum(axis=2)                 # |E_kj|_1
    def H(self, p):
        """floating-point H(p) (numerics only)"""
        return self.coef * np.exp(1j * (self.E @ np.asarray(p, dtype=float)))
    def C(self, p):
        if self.CE is None: return 0.0
        return np.exp(1j * (self.CE @ np.asarray(p, dtype=float))).sum()

def _circ_exp(e0, e1, e2):
    return np.array([[e0, e1, e2], [e2, e0, e1], [e1, e2, e0]])

def szollosi():
    z = np.zeros(4, dtype=int)
    a, b, c = z, np.array([-1, 0, 0, 0]), np.array([-1, -1, 0, 0])      # 1, conj x, conj(xy)
    d, e, f = z, np.array([0, 0, -1, 0]), np.array([0, 0, -1, -1])      # 1, conj u, conj(uv)
    EA = _circ_exp(a, b, c); EB = _circ_exp(d, e, f)
    EBs = _circ_exp(-d, -f, -e); EAs = _circ_exp(-a, -c, -b)            # B* = circ(d*, f*, e*), A* likewise
    E = np.zeros((6, 6, 4), dtype=int)
    E[:3, :3] = EA; E[:3, 3:] = EB; E[3:, :3] = EBs; E[3:, 3:] = EAs
    coef = np.ones((6, 6)); coef[3:, 3:] = -1.0
    CE = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [-1, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0, -1, -1]])
    return Family('szollosi', 4, E, coef, np.zeros((6, 6)), CE)

def fourier():
    """sqrt6 F(a,b), Jaming et al. (2009) eq. (3): rows
       [1, 1, 1, 1, 1, 1]
       [1, -w^2 x, w, -x, w^2, -w x]
       [1, w y, w^2, y, w, w^2 y]
       [1, -1, 1, -1, 1, -1]
       [1, w^2 x, w, x, w^2, w x]
       [1, -w y, w^2, -y, w, -w^2 y]      w = exp(2 pi i/3), x = e^{i s}, y = e^{i t}, p = (s, t)."""
    # coefficient as exponent of zeta = exp(i pi/3) (sixth root of unity): -1 = zeta^3, w = zeta^2, w^2 = zeta^4
    Z = np.array([[0, 0, 0, 0, 0, 0],
                  [0, 3 + 4, 2, 3, 4, 3 + 2],
                  [0, 2, 4, 0, 2, 4],
                  [0, 3, 0, 3, 0, 3],
                  [0, 4, 2, 0, 4, 2],
                  [0, 3 + 2, 4, 3, 2, 3 + 4]]) % 6
    X = np.array([[0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 1], [0, 0, 0, 0, 0, 0]])
    Y = np.array([[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 1],
                  [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 1, 0, 1, 0, 1]])
    E = np.stack([X, Y], axis=2)
    coef, err = sixth_root_table(Z)
    fam = Family('fourier', 2, E, coef, err, None)
    fam.Z = Z
    return fam

def sixth_root_table(Z):
    """float values of zeta^Z (zeta = e^{i pi/3}) and rigorous error bounds (mpmath interval arithmetic)"""
    import mpmath as mp
    mp.iv.dps = 30
    coef = np.empty(Z.shape, dtype=complex); err = np.empty(Z.shape)
    for idx, z in np.ndenumerate(Z):
        c = complex(np.cos(np.pi * z / 3), np.sin(np.pi * z / 3))
        if z == 0: c = 1 + 0j
        if z == 3: c = -1 + 0j
        th = mp.iv.pi * int(z) / 3
        er = abs(mp.iv.cos(th) - mp.iv.mpf(c.real)); ei = abs(mp.iv.sin(th) - mp.iv.mpf(c.imag))
        bound = float(er.b + ei.b)                           # |.|_2 <= |.|_1; float() rounds to nearest
        coef[idx] = c; err[idx] = float(np.nextafter(bound, np.inf))
    return coef, err
