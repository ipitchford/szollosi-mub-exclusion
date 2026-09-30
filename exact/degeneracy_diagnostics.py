"""
degeneracy_diagnostics.py -- NUMERICAL DIAGNOSTICS (not part of the proof) of the unbiased-vector sets at special
parameters of the Szollosi family, in 50-digit arithmetic (mpmath).  Referee comment 4.6: these observations are
reported as diagnostics; the certificate does not depend on them.

For a parameter point p (given through alpha and an ordering of the cubic roots) we
  1. locate solutions of |(H^* v)_j|^2 = 6 (j = 1..6), v = (1, e^{i theta_2..6}), by numpy multistart Newton;
  2. refine every solution with Newton / Gauss-Newton in 50-digit arithmetic and deduplicate at 1e-25;
  3. report the number of distinct refined solutions, the residuals, and the singular values of the 6x5 real
     Jacobian (a singular value ~1e-45 indicates a singular solution, which may be an isolated multiple root or
     lie on a positive-dimensional component -- the diagnostic does not decide which);
  4. along paths approaching the ray arg(alpha) = pi/6 and alpha = 0, report the minimum distance between distinct
     solutions (coalescence rate).
"""
import numpy as np, mpmath as mp, json, sys, time
mp.mp.dps = 60
rng = np.random.default_rng(2026)

def roots_mp(alpha):
    """roots of t^3 - alpha t^2 + conj(alpha) t - 1 (mpmath polyroots, 60 digits), sorted by argument in [0, 2pi)"""
    r = mp.polyroots([1, -alpha, mp.conj(alpha), -1], maxsteps=200, extraprec=200)
    return sorted(r, key=lambda t: float(mp.arg(t)) % (2 * np.pi))

def H_mp(alpha):
    x, y, z = roots_mp(alpha); u, v, w = roots_mp(-alpha)
    def circ(a0, a1, a2): return [[a0, a1, a2], [a2, a0, a1], [a1, a2, a0]]
    A = circ(1, mp.conj(x), z); B = circ(1, mp.conj(u), w)
    M = mp.matrix(6, 6)
    for i in range(3):
        for j in range(3):
            M[i, j] = A[i][j]; M[i, j + 3] = B[i][j]; M[i + 3, j] = mp.conj(B[j][i]); M[i + 3, j + 3] = -mp.conj(A[j][i])
    return M, (x, y, z, u, v, w)

def solve_np(Hn, nstart=40000):
    Hs = Hn.conj().T
    T = rng.uniform(0, 2 * np.pi, (nstart, 5))
    for it in range(60):
        V = np.exp(1j * np.concatenate([np.zeros((len(T), 1)), T], axis=1)); W = V @ Hs.T
        G = (np.abs(W) ** 2 - 6.0)                                      # 6 equations
        J = 2 * np.real(np.conj(W)[:, :, None] * Hs[None, :, 1:] * (1j * V[:, None, 1:]))   # 6 x 5
        with np.errstate(all='ignore'):
            dT = np.linalg.lstsq(J[0], G[0], rcond=None)[0] if False else np.einsum('nij,nj->ni', np.linalg.pinv(J), G)
        dT = np.nan_to_num(dT); st = np.clip(np.abs(dT).max(axis=1, keepdims=True), 1e-300, None)
        T = (T - np.minimum(1, 0.5 / st) * dT) % (2 * np.pi)
    V = np.exp(1j * np.concatenate([np.zeros((len(T), 1)), T], axis=1))
    res = np.abs(np.abs(V @ Hs.T) ** 2 - 6).max(axis=1)
    T = T[res < 1e-9]
    out = []
    for t in T:
        if not out or np.abs(np.exp(1j * t) - np.exp(1j * np.array(out))).max(axis=1).min() > 1e-6: out.append(t)
    return np.array(out)

def jac_mp(Hs, v, Wv):
    J = mp.matrix(6, 5)
    for j in range(6):
        for k in range(5):
            J[j, k] = 2 * mp.re(mp.conj(Wv[j]) * Hs[j, k + 1] * 1j * v[k + 1])
    return J

def refine_mp(H, t0, iters=400):
    """Gauss-Newton in 60-digit arithmetic (quadratic convergence at simple roots, linear at singular ones)"""
    th = [mp.mpf(float(a)) for a in t0]
    Hs = H.H
    for it in range(iters):
        v = mp.matrix([1] + [mp.expj(a) for a in th]); Wv = Hs * v
        G = mp.matrix([abs(Wv[j]) ** 2 - 6 for j in range(6)])
        J = jac_mp(Hs, v, Wv)
        JT = J.T
        lam = mp.norm(G) + mp.mpf(10) ** (-100)          # Levenberg-Marquardt damping (keeps singular roots solvable)
        try:
            d = mp.lu_solve(JT * J + mp.eye(5) * lam, JT * G)
        except ZeroDivisionError:
            break
        th = [th[k] - d[k] for k in range(5)]
        if mp.norm(d) < mp.mpf(10) ** (-56): break
    v = mp.matrix([1] + [mp.expj(a) for a in th]); Wv = Hs * v
    res = max(abs(abs(Wv[j]) ** 2 - 6) for j in range(6))
    sv = mp.svd_r(jac_mp(Hs, v, Wv), compute_uv=False)
    return th, res, sorted([sv[i] for i in range(len(sv))])

def solutions(alpha, nstart=40000, extra=None):
    H, roots = H_mp(alpha)
    Hn = np.array([[complex(H[i, j]) for j in range(6)] for i in range(6)])
    sols = solve_np(Hn, nstart)
    if extra is not None and len(extra): sols = np.concatenate([sols, extra])
    ref = []
    for s in sols:
        th, res, sv = refine_mp(H, s)
        if res > mp.mpf(10) ** (-40): continue
        vec = [mp.expj(a) for a in th]
        if any(max(abs(vec[k] - r[0][k]) for k in range(5)) < mp.mpf(10) ** (-18) for r in ref): continue
        ref.append((vec, res, sv, [float(t) % (2 * np.pi) for t in th]))
    return ref

def clusters(ref, thr):
    """single-linkage clusters of solutions at distance <= thr; returns sorted cluster sizes"""
    n = len(ref); parent = list(range(n))
    def f(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            if max(abs(ref[i][0][k] - ref[j][0][k]) for k in range(5)) <= thr: parent[f(i)] = f(j)
    from collections import Counter
    return sorted(Counter(f(i) for i in range(n)).values(), reverse=True)

def analyse(label, alpha, nstart=40000, extra=None, cluster_thr=None):
    t0 = time.time()
    ref = solutions(alpha, nstart, extra)
    smin = [r[2][0] for r in ref]
    tiny = mp.mpf(10) ** -20
    dmin = mp.inf
    for i in range(len(ref)):
        for j in range(i + 1, len(ref)):
            dmin = min(dmin, max(abs(ref[i][0][k] - ref[j][0][k]) for k in range(5)))
    nons = [s for s in smin if s >= tiny]
    out = {'label': label, 'alpha': [str(mp.nstr(mp.re(alpha), 20)), str(mp.nstr(mp.im(alpha), 20))],
           'distinct_solutions': len(ref), 'max_residual': float(max(r[1] for r in ref)) if ref else None,
           'singular_solutions': sum(1 for s in smin if s < tiny),
           'log10_min_smallest_singular_value_nonsingular': round(float(mp.log10(min(nons))), 3) if nons else None,
           'min_distance_between_solutions': float(dmin) if len(ref) > 1 else None,
           'time_s': round(time.time() - t0, 1)}
    if cluster_thr is not None: out['cluster_sizes_at_threshold'] = [float(cluster_thr), clusters(ref, cluster_thr)]
    print(json.dumps(out), flush=True)
    return out, ref

if __name__ == '__main__':
    results = []
    ray = mp.mpf(1) / 2 + 1j / (2 * mp.sqrt(3))                   # exact point of the ray arg = pi/6, |alpha| = 1/sqrt3
    results.append(analyse('generic alpha = 0.7 + 0.1 i', mp.mpc('0.7', '0.1'))[0])
    results.append(analyse('boundary point alpha = 1', mp.mpc(1, 0))[0])
    near = None
    for d in ['1e-2', '1e-4', '1e-6']:
        a = abs(ray) * mp.expj(mp.pi / 6 - mp.mpf(d))
        o, near = analyse(f'near the ray: arg alpha = pi/6 - {d}, |alpha| = 1/sqrt3', a, cluster_thr=10 * mp.mpf(d))
        results.append(o)
    extra = np.array([r[3] for r in near])
    results.append(analyse('ray point alpha = 1/2 + i/(2 sqrt3) (exact)', ray, extra=extra)[0])
    for r in ['1e-2', '1e-4', '1e-6']:
        o, near = analyse(f'near 0: alpha = {r} exp(i pi/12)', mp.mpf(r) * mp.expj(mp.pi / 12), cluster_thr=10 * mp.mpf(r))
        results.append(o)
    extra = np.array([r[3] for r in near])
    results.append(analyse('alpha = 0 (exact)', mp.mpc(0, 0), extra=extra)[0])
    json.dump(results, open('degeneracy_diagnostics.json', 'w'), indent=1)
