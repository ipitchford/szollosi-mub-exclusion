"""
driver.py -- full generator run over the selected parameter boxes of a family (GENERATOR).

Roots: selected boxes at level LROOT; children at level LCHILD selected by the same rigorous selector.
For each root: modulus-only root pass (phase levels 4..7), numerical solutions at the root, then each child is
certified with solve.child_box (first Lmax=9, then Lmax=10); a child that still fails is subdivided (selector
applied to the grandchildren) down to level LMAXP.

Record states (referee comment 4.5): every record has 'status' in {'certified', 'failed', 'refined'};
  certified : the box's certificate succeeded (its hulls are stored);
  failed    : the certificate failed at the deepest level (proof incomplete -- reported, never counted as done);
  refined   : the certificate failed and the box was subdivided; the box counts as covered only if every
              selected child is covered (checked by the checker's cover verification, never assumed here).
On resume, ONLY 'certified' records are skipped; 'refined' and 'failed' boxes are re-processed.

Usage: python3 driver.py FAMILY NPROC [k1,k2/m stride] [limit]      (results -> $MUB6_RESULTS)
"""
import numpy as np, json, os, sys, time
from multiprocessing import Pool
import family as fam_mod, core, solve, selector
from core import ParamBox

LROOT, LCHILD, LMAXP = 7, 9, 12
BASE = os.environ.get('MUB6_BASE', os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(BASE, 'results', os.environ.get('MUB6_RESULTS', 'results.jsonl'))
FAM = None

def get_family(name):
    return fam_mod.szollosi() if name == 'szollosi' else fam_mod.fourier()

def selector_test(fam):
    return selector.szollosi_test if fam.name == 'szollosi' else selector.fourier_test

def record(res):
    os.makedirs(os.path.dirname(RESULTS), exist_ok=True)
    with open(RESULTS, 'a') as f:
        f.write(json.dumps(res) + '\n')

def certified_set():
    s = set()
    if os.path.exists(RESULTS):
        for line in open(RESULTS):
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get('status') == 'certified' and r.get('ok') is True:
                s.add((r['lev'], tuple(r['idx'])))
    return s

def certify_box(P, idx7, L7, sols, done):
    """returns (n_certified_leaves, n_failed_leaves) of the subtree of P"""
    key = (P.lev, tuple(P.idx.tolist()))
    if key in done:
        return (1, 0)
    t0 = time.time()
    info = None
    for Lmax, Ld in [(9, 7), (10, 8)]:
        try:
            info = solve.child_box(P, idx7, L7, sols, Lmax=Lmax, Ld=Ld)
        except Exception as e:
            info = {'lev': P.lev, 'idx': P.idx.tolist(), 'ok': False, 'reason': 'exception: ' + repr(e)}
        sols = info.pop('sols', sols)
        info['Lmax'] = Lmax
        if info.get('ok'): break
    info['time'] = round(time.time() - t0, 1)
    if info.get('ok'):
        info['status'] = 'certified'; record(info); return (1, 0)
    if P.lev >= LMAXP:
        info['status'] = 'failed'; record(info); return (0, 1)
    info['status'] = 'refined'; info.pop('hull_lo', None); info.pop('hull_hi', None); record(info)
    ch = np.array([c.idx for c in P.children()])
    keep = selector_test(P.fam)(ch, P.lev + 1)
    nok = nf = 0
    for ci in ch[keep]:
        a, b = certify_box(ParamBox(P.fam, ci, P.lev + 1), idx7, L7, sols, done)
        nok += a; nf += b
    return (nok, nf)

def process_root(args):
    fam_name, root_idx, done = args
    fam = get_family(fam_name)
    P = ParamBox(fam, root_idx, LROOT)
    t0 = time.time()
    k = LCHILD - LROOT
    offs = np.indices((2 ** k,) * fam.m).reshape(fam.m, -1).T
    cidx = P.idx * 2 ** k + offs
    cidx = cidx[selector_test(fam)(cidx, LCHILD)]
    if len(cidx) == 0 or all((LCHILD, tuple(c.tolist())) in done for c in cidx):
        return {'root': root_idx.tolist(), 'children': int(len(cidx)), 'skipped': True, 'time': 0.0}
    idx7, L7, counts = solve.root_pass(P)
    th = solve.on_surface_point(P)
    sols = solve.multistart_solutions(fam.H(th), nstart=20000, seed=0)
    nok = nf = 0
    for ci in cidx:
        a, b = certify_box(ParamBox(fam, ci, LCHILD), idx7, L7, sols, done)
        nok += a; nf += b
    return {'root': root_idx.tolist(), 'children': int(len(cidx)), 'ok': nok, 'fail': nf, 'n7': int(len(idx7)),
            'root_counts': counts, 'nsol': int(len(sols)), 'time': round(time.time() - t0, 1)}

def roots_for(fam):
    roots, _ = selector.select(selector_test(fam), fam.m, LROOT)
    return roots

if __name__ == "__main__":
    fam_name = sys.argv[1]; nproc = int(sys.argv[2])
    stride = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != '-' else None
    limit = int(sys.argv[4]) if len(sys.argv) > 4 else None
    fam = get_family(fam_name)
    roots = roots_for(fam)
    TWO_PI = 2 * np.pi
    def key(r):
        th = (r + 0.5) * TWO_PI / 2 ** LROOT
        if fam.name == 'szollosi':
            x, y = np.exp(1j * th[0]), np.exp(1j * th[1]); return -abs(x + y + np.conj(x * y))
        return tuple(r)
    roots = sorted(roots, key=key)
    if stride:
        ks, mm = stride.split('/'); mm = int(mm); kk = set(int(x) for x in ks.split(','))
        roots = [r for i, r in enumerate(roots) if i % mm in kk]
    if os.environ.get('MUB6_REVERSE') == '1': roots = roots[::-1]      # tail-first (load balancing between machines)
    if limit: roots = roots[:limit]
    print(f"{fam_name}: {len(roots)} roots at level {LROOT}; {nproc} processes", flush=True)
    core.table17()
    done = certified_set()
    print(f"resume: {len(done)} certified records found", flush=True)
    t0 = time.time()
    with Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(process_root, [(fam_name, r, done) for r in roots])):
            print(f"[{i+1}/{len(roots)}] {r}  elapsed {time.time()-t0:.0f}s", flush=True)
    print("done", flush=True)
