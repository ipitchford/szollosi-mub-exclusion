"""
core.py -- rigorous enclosures, exclusion tests (modulus, mean-value modulus, parametric Krawczyk) for a general
admissible family (family.py).  GENERATOR side.  Every radius is passed through rnd.up (Lemma R); every exclusion
is a strict comparison lower > upper of rigorously rounded quantities (rnd.nd / rnd.nu).

Grids: centres of phase and parameter boxes lie on G17 = (2 pi / 2^17) Z; box end-points on G16.
A phase box is (lo, hi) in Z^5 on G16 with lo <= hi < lo + 2^16; its centre index on G17 is lo + hi and its
angular half-width is (hi - lo) pi / 2^16.
"""
import numpy as np, os
from rnd import U, TAB_ERR, up, nu, nd, SQ6, angle_up, PI_UP

GL = 17; NG = 2 ** GL; NG16 = 2 ** 16
_T17 = None

def table17():
    """unit table T[m] ~ exp(2 pi i m / 2^17); loaded from disk and verified with Arb on every load"""
    global _T17
    if _T17 is not None: return _T17
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'table17.npy')
    t = np.load(fn)
    if os.environ.get('MUB6_SKIP_TABLE_VERIFY') != '1':
        _verify17(t)
    _T17 = t
    return t

def _verify17(t):
    import flint
    flint.ctx.prec = 120
    pi = flint.arb.pi(); bound = flint.arb(2.0 ** -53); bad = 0
    for m in range(NG):
        th = 2 * pi * flint.arb(m) / flint.arb(NG)
        er = (th.cos() - flint.arb(float(t[m].real))).abs_upper()
        ei = (th.sin() - flint.arb(float(t[m].imag))).abs_upper()
        if not (bool(er <= bound) and bool(ei <= bound)): bad += 1
    if bad: raise RuntimeError(f"table17: {bad} entries exceed the rounding bound")

def unit_at(cidx):
    return table17()[np.asarray(cidx) % NG]

# ------------------------------------------------------------------------------------------------ disc primitives
def dmul(am, ar, bm, br):
    """product of complex discs (am, ar) (bm, br): radius |am| br + ar |bm| + ar br + 4U |am||bm| (complex product
    rounding <= sqrt5 u |am||bm|; abs() within 2u), made rigorous by up()."""
    m = am * bm
    aa = np.abs(am); ab = np.abs(bm)
    r = up(aa * br + ar * ab + ar * br + 4 * U * aa * ab)
    return m, r

def dabs(m, r):
    """interval (a, ar) containing |z| for z in disc (m, r): a = fl|m| (hypot, <= 2u relative error)"""
    a = np.abs(m)
    return a, up(r + 4 * U * a)

def imatmul_point(Y, Bm, Br):
    """(point matrix Y) x (interval matrix Bm +- Br), batched: midpoint Y Bm, radius |Y| Br + (p+1) U |Y||Bm|"""
    Cm = Y @ Bm
    p = Y.shape[-1]
    Ya = np.abs(Y)
    Cr = up(Ya @ Br + (p + 1) * U * (Ya @ np.abs(Bm)))
    return Cm, Cr

# ------------------------------------------------------------------------------------------------ parameter boxes
class ParamBox:
    """dyadic parameter box of level lev (<= 16) and index idx (m,) for family fam.
    Angles p_i in [idx_i, idx_i + 1] * 2 pi / 2^lev; centre index on G17: (2 idx + 1) 2^(16 - lev)."""
    def __init__(self, fam, idx, lev):
        assert lev <= 16
        self.fam = fam; self.idx = np.asarray(idx, dtype=np.int64); self.lev = lev
        self.half_g17 = 2 ** (16 - lev)                           # half-width in G17 units
        self.hw = float(angle_up(self.half_g17, 17))             # >= pi / 2^lev
        self.cidx = (2 * self.idx + 1) * self.half_g17
        E = fam.E
        ph = (E @ self.cidx) % NG                                  # <E_kj, p_c> in G17 units (exact integers)
        T = unit_at(ph)
        Hm = fam.coef * T
        # |coef T - coef_exact T_exact| <= |coef| TAB_ERR + err |T| + err TAB_ERR ;  rounding of the product
        Hr0 = up(np.abs(fam.coef) * TAB_ERR + fam.coef_err * np.abs(T) + fam.coef_err * TAB_ERR
                 + 4 * U * np.abs(fam.coef) * np.abs(T))
        self.Hm, self.Hr0 = Hm, Hr0
        self.Hc = np.conj(Hm); self.Habs = np.abs(Hm)
        if fam.CE is not None:
            ct = (fam.CE @ self.cidx) % NG
            terms = unit_at(ct)
            self.Cterms = terms
            self.S_c = terms.sum()
            self.S_cr = float(up(len(terms) * TAB_ERR + len(terms) * U * np.abs(terms).sum()))
            self.D_c = np.array([1j * (fam.CE[:, i] * terms).sum() for i in range(fam.m)])
            best, bestdet = None, -1.0
            for a in range(fam.m):
                for b in range(a + 1, fam.m):
                    M = np.array([[self.D_c[a].real, self.D_c[b].real], [self.D_c[a].imag, self.D_c[b].imag]])
                    d = abs(np.linalg.det(M))
                    if d > bestdet: bestdet, best = d, (a, b)
            self.dep = best
        else:
            self.dep = ()
        self.free = tuple(i for i in range(fam.m) if i not in self.dep)
    def children(self):
        offs = np.indices((2,) * self.fam.m).reshape(self.fam.m, -1).T
        return [ParamBox(self.fam, self.idx * 2 + o, self.lev + 1) for o in offs]

# ------------------------------------------------------------------------------------------------ phase boxes
def box_geometry(lo, hi):
    """centres V (N,6) (V[:,0] = 1) and radii rho (N,6): rho_k >= |exp(i theta_k) - V_k| for theta in the box
    (chord <= arc <= angular half-width, plus table error); rho[:,0] = 0."""
    lo = np.asarray(lo); hi = np.asarray(hi)
    c = (lo + hi) % NG
    V = np.ones((lo.shape[0], 6), dtype=complex); V[:, 1:] = unit_at(c)
    rho = np.zeros((lo.shape[0], 6))
    rho[:, 1:] = up(angle_up(hi - lo, 17) + TAB_ERR)          # (hi-lo) G17 units = (hi-lo) pi / 2^16
    return V, rho

def half_widths(lo, hi):
    """angular half-widths (upper bounds) of the phase coordinates 2..6"""
    return angle_up(np.asarray(hi) - np.asarray(lo), 17)

def dyadic_boxes(idx, L):
    s = 2 ** (16 - L)
    return idx * s, (idx + 1) * s

# ------------------------------------------------------------------------------------------------ enclosures
def _delta(P, hwv, weight=None):
    """|H_kj(p) - Hm_kj| <= sum_i |E_kji| hwv_i + Hr0_kj  (p within hwv of the centre); weight multiplies"""
    E = np.abs(P.fam.E).astype(float)
    return up(E @ hwv + P.Hr0)

def S_disc(V, rho, P, hwv):
    delta = _delta(P, hwv)
    Sm = V @ P.Hc
    Sr = up(rho @ P.Habs + delta.sum(axis=0)[None, :] + 11 * U * (np.abs(V) @ P.Habs))
    return Sm, Sr, delta

def _modulus_gap_positive(a, ar):
    """rigorous test | a_exact - sqrt6 | > ar  where a = fl|.| (see rnd.py)"""
    d = a - SQ6
    return nd(np.abs(d) * (1 - 2 * U)) > up(ar + U * SQ6)

def modulus_exclude(V, rho, P):
    hwv = np.full(P.fam.m, P.hw)
    Sm, Sr, _ = S_disc(V, rho, P, hwv)
    a, ar = dabs(Sm, Sr)
    return _modulus_gap_positive(a, ar).any(axis=1)

def _dS(V, rho, P, i, delta):
    """disc of dS_j/dp_i = -i sum_k E_kji conj(H_kj) v_k"""
    Ei = P.fam.E[:, :, i].astype(float)
    m = -1j * (V @ (Ei * P.Hc))
    r = up(rho @ (np.abs(Ei) * P.Habs) + (np.abs(Ei) * delta).sum(axis=0)[None, :]
           + 11 * U * (np.abs(V) @ (np.abs(Ei) * P.Habs)))
    return m, r

def _d2S(V, rho, P, a, b, delta):
    Eab = (P.fam.E[:, :, a] * P.fam.E[:, :, b]).astype(float)
    m = -(V @ (Eab * P.Hc))
    r = up(rho @ (np.abs(Eab) * P.Habs) + (np.abs(Eab) * delta).sum(axis=0)[None, :]
           + 11 * U * (np.abs(V) @ (np.abs(Eab) * P.Habs)))
    return m, r

def _Cderiv(P, hwv, idxs):
    """disc of the derivative of C w.r.t. the parameter angles idxs, over the parameter ranges hwv"""
    CE = P.fam.CE
    coef = np.ones(CE.shape[0])
    for i in idxs: coef = coef * CE[:, i]
    fac = (1j) ** len(idxs)
    var = up(np.abs(CE) @ hwv)                                 # |term(p) - term(c)| <= |e_t|_1-weighted hw
    m = fac * (coef * P.Cterms).sum()
    r = float(up((np.abs(coef) * (TAB_ERR + var)).sum() + 6 * U * (np.abs(coef) * np.abs(P.Cterms)).sum()))
    return m, r

def mv_modulus_exclude(V, rho, P):
    """mean-value modulus test (see paper): | |S_j(x)| - |S_j(c)| | <= sum_k G_jk rho_k + sum_i G_ji hw"""
    m_ = P.fam.m
    N = V.shape[0]
    hwv = np.full(m_, P.hw)
    Sm, Sr, delta = S_disc(V, rho, P, hwv)
    a = np.abs(Sm)
    lower = nd(nd(a * (1 - 4 * U)) - Sr)                        # rigorous lower bound of |S_j| on the box
    pos = lower > 0
    rho0 = np.full_like(rho, TAB_ERR); rho0[:, 0] = 0.0
    Sm0, Sr0, _ = S_disc(V, rho0, P, np.zeros(m_))
    a0, ar0 = dabs(Sm0, Sr0)
    safe_lower = np.where(pos, lower, 1.0)
    var = np.zeros((N, 6))
    for kk in range(1, 6):
        t1m, t1r = dmul(np.conj(Sm), Sr, P.Hc[kk][None, :], delta[kk][None, :])
        t2m, t2r = dmul(t1m, t1r, (1j * V[:, kk])[:, None], rho[:, kk][:, None])
        num = up(np.abs(t2m.real) + t2r)
        G = np.where(pos, nu(num / safe_lower), 1.0)
        var = up(var + G * rho[:, kk][:, None])
    for i in range(m_):
        dSm, dSr = _dS(V, rho, P, i, delta)
        t1m, t1r = dmul(np.conj(Sm), Sr, dSm, dSr)
        num = up(np.abs(t1m.real) + t1r)
        G = np.where(pos, nu(num / safe_lower), np.abs(P.fam.E[:, :, i]).sum(axis=0)[None, :].astype(float))
        var = up(var + G * P.hw)
    return _modulus_gap_positive(a0, up(ar0 + var)).any(axis=1)

# ------------------------------------------------------------------------------------------------ Krawczyk
def krawczyk(V, rho, hw5, P, return_K=False):
    """Parametric augmented Krawczyk exclusion test (second-order Taylor in the free parameters, paper Sec. 6.3).
    Unknowns w = (theta_2..theta_6, dependent parameter angles); parameters q = free angles.
    hw5: (N,5) angular half-widths of the phase box (upper bounds).  Returns boolean exclusion mask."""
    fam = P.fam; m_ = fam.m
    N = V.shape[0]
    hw = P.hw; dep, free = P.dep, P.free
    nd_ = len(dep); n = 5 + nd_
    hwP = np.full(m_, hw); hwQ = np.zeros(m_); hwQ[list(free)] = hw; hw0 = np.zeros(m_)
    hwD = np.zeros(m_); hwD[list(dep)] = hw
    rho0 = np.full_like(rho, TAB_ERR); rho0[:, 0] = 0.0
    # ---------- value at (w~; c)
    Sm, Src, delta0 = S_disc(V, rho0, P, hw0)
    a, ar = dabs(Sm, Src)
    Fm_all = a * a - 6.0
    Fr_all = up(2 * a * ar + ar * ar + 4 * U * a * a + 2 * U * np.abs(Fm_all))
    F0m = np.zeros((N, n)); F0r = np.zeros((N, n))
    F0m[:, :5] = Fm_all[:, 1:]; F0r[:, :5] = Fr_all[:, 1:]
    if nd_:
        F0m[:, 5] = P.S_c.real; F0m[:, 6] = P.S_c.imag; F0r[:, 5] = P.S_cr; F0r[:, 6] = P.S_cr
    nq = len(free)
    # ---------- first derivatives in q at (w~; c)
    Fq_m = np.zeros((N, n, nq)); Fq_r = np.zeros((N, n, nq))
    for t, i in enumerate(free):
        dSm, dSr = _dS(V, rho0, P, i, delta0)
        pm, pr = dmul(np.conj(Sm), Src, dSm, dSr)
        Fq_m[:, :5, t] = 2 * pm.real[:, 1:]; Fq_r[:, :5, t] = up(2 * pr[:, 1:] + 4 * U * np.abs(pm.real[:, 1:]))
        if nd_:
            cm, cr = _Cderiv(P, hw0, (i,))
            Fq_m[:, 5, t] = cm.real; Fq_m[:, 6, t] = cm.imag; Fq_r[:, 5, t] = cr; Fq_r[:, 6, t] = cr
    # ---------- second derivatives in q over Q (theta = centre)
    SmQ, SrQ, deltaQ = S_disc(V, rho0, P, hwQ)
    Fqq_m = np.zeros((N, n, nq, nq)); Fqq_r = np.zeros((N, n, nq, nq))
    dSQ = {i: _dS(V, rho0, P, i, deltaQ) for i in free}
    for ta, ia in enumerate(free):
        for tb, ib in enumerate(free):
            d2m, d2r = _d2S(V, rho0, P, ia, ib, deltaQ)
            p1m, p1r = dmul(np.conj(dSQ[ia][0]), dSQ[ia][1], dSQ[ib][0], dSQ[ib][1])
            p2m, p2r = dmul(np.conj(SmQ), SrQ, d2m, d2r)
            val_m = 2 * (p1m.real + p2m.real)
            val_r = up(2 * (p1r + p2r) + 4 * U * (np.abs(p1m.real) + np.abs(p2m.real)))
            Fqq_m[:, :5, ta, tb] = val_m[:, 1:]; Fqq_r[:, :5, ta, tb] = val_r[:, 1:]
            if nd_:
                cm, cr = _Cderiv(P, hwQ, (ia, ib))
                Fqq_m[:, 5, ta, tb] = cm.real; Fqq_m[:, 6, ta, tb] = cm.imag
                Fqq_r[:, 5, ta, tb] = cr; Fqq_r[:, 6, ta, tb] = cr
    # ---------- Jacobian over X x D (free params at centre)
    SmX, SrX, deltaD = S_disc(V, rho, P, hwD)
    Jm = np.zeros((N, n, n)); Jr = np.zeros((N, n, n))
    for kk in range(1, 6):
        t1m, t1r = dmul(np.conj(SmX[:, 1:]), SrX[:, 1:], P.Hc[kk, 1:][None, :], deltaD[kk, 1:][None, :])
        t2m, t2r = dmul(t1m, t1r, (1j * V[:, kk])[:, None], rho[:, kk][:, None])
        Jm[:, :5, kk - 1] = 2 * t2m.real; Jr[:, :5, kk - 1] = up(2 * t2r + 4 * U * np.abs(t2m.real))
    for t, i in enumerate(dep):
        dSm, dSr = _dS(V, rho, P, i, deltaD)
        pm, pr = dmul(np.conj(SmX), SrX, dSm, dSr)
        Jm[:, :5, 5 + t] = 2 * pm.real[:, 1:]; Jr[:, :5, 5 + t] = up(2 * pr[:, 1:] + 4 * U * np.abs(pm.real[:, 1:]))
        cm, cr = _Cderiv(P, hwD, (i,))
        Jm[:, 5, 5 + t] = cm.real; Jm[:, 6, 5 + t] = cm.imag; Jr[:, 5, 5 + t] = cr; Jr[:, 6, 5 + t] = cr
    # ---------- q-derivatives of the Jacobian over X x P
    SmP, SrP, deltaP = S_disc(V, rho, P, hwP)
    dSP = {i: _dS(V, rho, P, i, deltaP) for i in range(m_)}
    Jq_m = np.zeros((N, n, n, nq)); Jq_r = np.zeros((N, n, n, nq))
    for t, i in enumerate(free):
        Ei = fam.E[:, :, i].astype(float)
        for kk in range(1, 6):
            dHc = (-1j * Ei[kk, 1:] * P.Hc[kk, 1:])[None, :]
            dHr = up(np.abs(Ei[kk, 1:]) * deltaP[kk, 1:])[None, :]
            a1m, a1r = dmul(np.conj(dSP[i][0][:, 1:]), dSP[i][1][:, 1:], P.Hc[kk, 1:][None, :], deltaP[kk, 1:][None, :])
            a2m, a2r = dmul(np.conj(SmP[:, 1:]), SrP[:, 1:], dHc, dHr)
            sm = a1m + a2m; sr = up(a1r + a2r + 2 * U * np.abs(sm))
            bm, br = dmul(sm, sr, (1j * V[:, kk])[:, None], rho[:, kk][:, None])
            Jq_m[:, :5, kk - 1, t] = 2 * bm.real; Jq_r[:, :5, kk - 1, t] = up(2 * br + 4 * U * np.abs(bm.real))
        for s_, d_ in enumerate(dep):
            d2m, d2r = _d2S(V, rho, P, i, d_, deltaP)
            a1m, a1r = dmul(np.conj(dSP[i][0]), dSP[i][1], dSP[d_][0], dSP[d_][1])
            a2m, a2r = dmul(np.conj(SmP), SrP, d2m, d2r)
            sm = a1m + a2m; sr = up(a1r + a2r + 2 * U * np.abs(sm))
            Jq_m[:, :5, 5 + s_, t] = 2 * sm.real[:, 1:]; Jq_r[:, :5, 5 + s_, t] = up(2 * sr[:, 1:] + 4 * U * np.abs(sm.real[:, 1:]))
            cm, cr = _Cderiv(P, hwP, (i, d_))
            Jq_m[:, 5, 5 + s_, t] = cm.real; Jq_m[:, 6, 5 + s_, t] = cm.imag
            Jq_r[:, 5, 5 + s_, t] = cr; Jq_r[:, 6, 5 + s_, t] = cr
    # ---------- point Jacobian at the centre and Y (numerics only; any Y is admissible)
    Jc = np.zeros((N, n, n))
    for kk in range(1, 6):
        Jc[:, :5, kk - 1] = 2 * np.real(np.conj(Sm[:, 1:]) * P.Hc[kk, 1:][None, :] * (1j * V[:, kk])[:, None])
    for t, i in enumerate(dep):
        dSc = -1j * (V @ (fam.E[:, :, i].astype(float) * P.Hc))
        Jc[:, :5, 5 + t] = 2 * np.real(np.conj(Sm) * dSc)[:, 1:]
        Jc[:, 5, 5 + t] = P.D_c[i].real; Jc[:, 6, 5 + t] = P.D_c[i].imag
    ok = np.ones(N, dtype=bool); Y = np.zeros((N, n, n))
    with np.errstate(all='ignore'):
        try:
            Y = np.linalg.inv(Jc)
        except np.linalg.LinAlgError:
            for k in range(N):
                try: Y[k] = np.linalg.inv(Jc[k])
                except np.linalg.LinAlgError: ok[k] = False
    ok &= np.isfinite(Y).all(axis=(1, 2)); Y[~ok] = 0.0
    # ---------- Y F~(w~; Q) = Y F0 + (Y Fq) dq + 1/2 dq^T (Y Fqq) dq
    YFm, YFr = imatmul_point(Y, F0m[:, :, None], F0r[:, :, None]); YFm, YFr = YFm[:, :, 0], YFr[:, :, 0]
    if nq:
        YQm, YQr = imatmul_point(Y, Fq_m, Fq_r)
        YQQm, YQQr = imatmul_point(Y, Fqq_m.reshape(N, n, nq * nq), Fqq_r.reshape(N, n, nq * nq))
        YFr = up(YFr + up(np.abs(YQm) + YQr).sum(axis=2) * hw + 0.5 * up(np.abs(YQQm) + YQQr).sum(axis=2) * (hw * hw))
    # ---------- E = I - Y J~(W;Q) = I - Y J(W;c) - (Y J_q) dq
    YJm, YJr = imatmul_point(Y, Jm, Jr)
    Em = np.eye(n)[None] - YJm
    Er = up(YJr + 2 * U * np.abs(Em))
    if nq:
        YJq_m, YJq_r = imatmul_point(Y, Jq_m.reshape(N, n, n * nq), Jq_r.reshape(N, n, n * nq))
        YJq_m = YJq_m.reshape(N, n, n, nq); YJq_r = YJq_r.reshape(N, n, n, nq)
        Er = up(Er + up(np.abs(YJq_m) + YJq_r).sum(axis=3) * hw)
    radX = np.concatenate([hw5, np.full((N, nd_), hw)], axis=1)          # half-widths of W
    EXr = up((up(np.abs(Em) + Er) @ radX[:, :, None])[:, :, 0])
    Krad = up(YFr + EXr)
    Kshift = -YFm
    excl = (np.abs(Kshift) > up(Krad + radX)).any(axis=1) & ok
    if return_K:
        return excl, Kshift, Krad, ok
    return excl
