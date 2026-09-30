"""trial run of the repaired generator on sample boxes: python3 tests/trial_gen.py boxes.json nproc"""
import sys, os, json, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gen'))
import family, core, solve, selector
def run(idx):
    fam = family.szollosi(); idx = np.array(idx)
    if not selector.szollosi_test(idx[None, :], 9)[0]:
        return {'idx': idx.tolist(), 'selected': False}
    t0 = time.time()
    root = core.ParamBox(fam, idx // 4, 7)
    idx7, L7, counts = solve.root_pass(root)
    t1 = time.time()
    sols = solve.multistart_solutions(fam.H(solve.on_surface_point(root)), nstart=20000, seed=0)
    P = core.ParamBox(fam, idx, 9)
    for Lmax, Ld in [(9, 7), (10, 8)]:
        info = solve.child_box(P, idx7, L7, sols, Lmax=Lmax, Ld=Ld); info.pop('sols', None)
        if info.get('ok'): break
    info.pop('hull_lo', None); info.pop('hull_hi', None)
    return {'idx': idx.tolist(), 'ok': info.get('ok'), 'reason': info.get('reason'), 'hulls': info.get('hulls'),
            'max_hw': info.get('max_hw'), 'cliques': info.get('cert', {}).get('cliques'), 'nsol': info.get('nsol'),
            'Lmax': Lmax, 'root_time': round(t1 - t0, 1), 'child_time': info.get('time'), 'n7': len(idx7), 'counts': counts}
if __name__ == '__main__':
    boxes = json.load(open(sys.argv[1]))
    core.table17()
    with Pool(int(sys.argv[2])) as pool:
        for r in pool.imap_unordered(run, boxes):
            print(json.dumps(r), flush=True)
