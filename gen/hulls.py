"""
hulls.py -- cyclic hulls and clustering of phase boxes (generator side).

cyclic_hull(lo, hi, N): per coordinate, the complement of the LARGEST gap of the union of the cyclic intervals
[lo_i, hi_i] (mod N).  Repair of referee comment 4.4: every interval is first reduced to [a, a+len] with
0 <= a < N and, if it wraps (a + len >= N), split into [a, N) and [0, a + len - N]; the pieces are sorted and
merged on [0, N); gaps are then read off exactly, including the wrap-around gap.  Contract (tested against a
brute-force check in tests/test_hulls.py): every point of every input interval lies in the returned interval.
"""
import numpy as np

NG16 = 2 ** 16

def cyclic_hull_1d(lo, hi, N=NG16):
    """lo, hi integer arrays (hi >= lo).  Returns (s, e) with e - s < N or (0, N): the smallest-length cyclic
    interval [s, e] (mod N) containing all [lo_i, hi_i]."""
    lo = np.asarray(lo, dtype=np.int64); hi = np.asarray(hi, dtype=np.int64)
    length = hi - lo
    if (length >= N - 1).any():
        return 0, N
    a = lo % N; b = a + length                         # closed interval [a, b], 0 <= a < N, b < a + N - 1
    pieces = []
    for ai, bi in zip(a.tolist(), b.tolist()):
        if bi < N:
            pieces.append((ai, bi))
        else:
            pieces.append((ai, N - 1)); pieces.append((0, bi - N))
    pieces.sort()
    merged = []
    for s, e in pieces:
        if merged and s <= merged[-1][1] + 1:          # integer points: adjacent intervals merge
            if e > merged[-1][1]: merged[-1][1] = e
        else:
            merged.append([s, e])
    if len(merged) == 1 and merged[0][0] == 0 and merged[0][1] == N - 1:
        return 0, N
    # gaps between consecutive merged pieces (points strictly between), and the wrap-around gap
    best_gap, best_after = -1, None
    for i in range(len(merged)):
        e = merged[i][1]
        s_next = merged[(i + 1) % len(merged)][0] + (N if i == len(merged) - 1 else 0)
        gap = s_next - e - 1                            # number of integer points not covered
        if gap > best_gap:
            best_gap, best_after = gap, i
    if best_gap <= 0:
        return 0, N
    # hull starts at the piece after the largest gap and ends at the piece before it
    start = merged[(best_after + 1) % len(merged)][0]
    end_piece = merged[best_after][1]
    end = end_piece if end_piece >= start else end_piece + N
    return start, end

def cyclic_hull(lo, hi, N=NG16):
    """per-coordinate cyclic hull of boxes (lo, hi) of shape (n, d)"""
    lo = np.atleast_2d(lo); hi = np.atleast_2d(hi)
    d = lo.shape[1]
    out_lo = np.zeros(d, dtype=np.int64); out_hi = np.zeros(d, dtype=np.int64)
    for k in range(d):
        out_lo[k], out_hi[k] = cyclic_hull_1d(lo[:, k], hi[:, k], N)
    return out_lo, out_hi

def contains(hlo, hhi, lo, hi, N=NG16):
    """boolean (n_boxes,): box [lo,hi] contained in cyclic hull [hlo,hhi] (per coordinate)"""
    L = hhi - hlo
    full = L >= N
    off = (lo - hlo) % N
    ok = full | ((off + (hi - lo)) <= L)
    return ok.all(axis=-1)

_NB5 = None
def components(N, rows, cols):
    """connected components of the undirected graph on N nodes with edges (rows[i], cols[i]): hook-and-compress
    label propagation in NumPy (no SciPy dependency).  Returns (number of components, labels 0..ncomp-1)."""
    lab = np.arange(N)
    rows = np.asarray(rows, dtype=np.int64); cols = np.asarray(cols, dtype=np.int64)
    while True:
        new = lab.copy()
        if len(rows):
            m = np.minimum(lab[rows], lab[cols])
            np.minimum.at(new, rows, m); np.minimum.at(new, cols, m)
        while True:                                   # pointer jumping
            nn = new[new]
            if np.array_equal(nn, new): break
            new = nn
        if np.array_equal(new, lab): break
        lab = new
    uniq, lab2 = np.unique(lab, return_inverse=True)
    return len(uniq), lab2.reshape(-1)

def cluster(idx, L):
    """connected components of touching dyadic boxes (same level L) (numerics only: the grouping never affects
    soundness, every surviving box ends up in some hull)"""
    global _NB5
    d = idx.shape[1]
    if _NB5 is None or _NB5.shape[1] != d:
        _NB5 = np.array([o for o in np.indices((3,) * d).reshape(d, -1).T - 1 if not (o == 0).all()])
    n = 2 ** L; N = len(idx)
    if N == 0: return []
    w = n ** np.arange(d - 1, -1, -1, dtype=np.int64)
    keys = idx.astype(np.int64) @ w
    order = np.argsort(keys); skeys = keys[order]
    rows = []; cols = []
    same = np.nonzero(skeys[1:] == skeys[:-1])[0]
    rows.append(order[same]); cols.append(order[same + 1])
    for o in _NB5:
        if not (tuple(o) > (0,) * d): continue
        nb = (idx + o) % n
        nkeys = nb.astype(np.int64) @ w
        pos = np.searchsorted(skeys, nkeys); pos[pos >= N] = 0
        hit = skeys[pos] == nkeys
        rows.append(np.nonzero(hit)[0]); cols.append(order[pos[hit]])
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    ncomp, lab = components(N, rows, cols)
    order2 = np.argsort(lab); lab_s = lab[order2]
    bounds = np.searchsorted(lab_s, np.arange(ncomp + 1))
    return [order2[bounds[c]:bounds[c + 1]] for c in range(ncomp)]
