"""
run_check.py -- independent replay of generator certificates (CHECKER driver).

Usage: python3 run_check.py FAMILY TABLE.npy OUT.jsonl NPROC GEN1.jsonl [GEN2.jsonl ...] [--roots k/m] [--rootlist roots.json] [--limit n]

1. The unit table is verified with mpmath interval arithmetic (ck.verify_table) and its SHA-256 recorded.
2. The cover tree is walked from level 3 with the CHECKER's selector: a kept box is 'covered' if it carries a
   certified generator record that the checker verifies, or if every kept child is covered (levels <= 12).
   The boxes that need verification are collected, grouped by their level-7 ancestor ('root').
3. For each root the checker runs its own modulus-only root pass (phase levels 4..7), then for each box
   verifies coverage (ck.verify_coverage) and the combinatorial certificate (ck.verify_graphs).
4. The cover tree is evaluated with the verification results; the summary states whether every point of the
   fundamental domain lies in a verified box.
No generator code is imported; generator records are read as plain JSON data.
"""
import sys, os, json, time, hashlib, platform, numpy as np
from multiprocessing import get_context
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ck
import receipt_binding as rb

LROOT, LTOP, LMAX = 7, 3, 12
G = {}

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def load_records(paths):
    recs = {}
    n = 0
    for p in paths:
        for line in open(p):
            try: r = json.loads(line)
            except Exception: continue
            n += 1
            if r.get('status') == 'certified' and r.get('ok') is True and 'hull_lo' in r:
                recs[(int(r['lev']), tuple(int(t) for t in r['idx']))] = (r['hull_lo'], r['hull_hi'])
    return recs, n

def walk(fam, recs):
    """boxes reached by the cover tree that carry a certified record; and leaves without records"""
    sel = ck.selector(fam); m = fam.m
    offs = np.indices((2,) * m).reshape(m, -1).T
    cur = np.indices((2 ** LTOP,) * m).reshape(m, -1).T.astype(np.int64)
    cur = cur[sel(cur, LTOP)]
    need, missing, kept_counts = [], [], {LTOP: len(cur)}
    lev = LTOP
    while len(cur) and lev <= LMAX:
        nxt = []
        for b in cur:
            key = (lev, tuple(int(t) for t in b))
            if key in recs: need.append(key)
            elif lev == LMAX: missing.append(key)
            else: nxt.append(b)
        if lev == LMAX or not nxt: break
        nxt = np.array(nxt, dtype=np.int64)
        ch = (nxt[:, None, :] * 2 + offs[None]).reshape(-1, m)
        cur = ch[sel(ch, lev + 1)]; lev += 1; kept_counts[lev] = len(cur)
    return need, missing, kept_counts

def evaluate(fam, verified):
    """cover evaluation: top-down from the kept level-3 boxes; a box is removed when it is verified, otherwise
    it is replaced by its kept children (checker selector); complete iff nothing remains after level LMAX"""
    ok_keys = {k for k, v in verified.items() if v is True}
    _, leftovers, _ = walk(fam, ok_keys)
    anc9 = sorted({tuple(int(t) >> (lev - 9) for t in idx) for (lev, idx) in leftovers})
    return len(leftovers) == 0, anc9

def init_worker(fam_name, table_path, table_sha, replay_context):
    if sha256(table_path) != table_sha:
        raise RuntimeError('table changed since verification')
    ck.load_table(table_path, verify=False)
    G['fam'] = ck.Fam(fam_name); G['fam'].install()
    G['context'] = replay_context

def process_root(args):
    root, items = args
    fam = G['fam']; t0 = time.time()
    pc, ph = ck.pbox(np.array(root), LROOT)
    levels, counts = ck.root_pass(pc, ph)
    out = []
    for (lev, idx, hlo, hhi) in items:
        if tuple(int(t) >> (lev - LROOT) for t in idx) != tuple(root):
            raise ValueError('box not below its root')
        t1 = time.time()
        hlo = np.asarray(hlo, dtype=np.int64).reshape(-1, 5); hhi = np.asarray(hhi, dtype=np.int64).reshape(-1, 5)
        ok_cov, st = ck.verify_coverage(fam, np.array(idx), lev, hlo, hhi, levels, L7=LROOT)
        ok_g, gi = ck.verify_graphs(hlo, hhi) if ok_cov else (False, {})
        hsha = hashlib.sha256(json.dumps([hlo.tolist(), hhi.tolist()]).encode()).hexdigest()[:16]
        out.append({'lev': lev, 'idx': list(idx), 'root': list(root), 'hulls_sha': hsha, 'ok': bool(ok_cov and ok_g), 'coverage': ok_cov,
                    'graphs': ok_g, 'nh': int(len(hlo)), 'cov_stats': st, 'graph_info': gi, 'time': round(time.time() - t1, 2),
                    'input_sha256': rb.binding(G['context'], (lev, idx), [hlo.tolist(), hhi.tolist()])})
    return {'root': list(root), 'root_counts': counts, 'boxes': out, 'time': round(time.time() - t0, 1), 'replay_context': G['context']}

def main():
    args = [a for a in sys.argv[1:]]
    clean = '--clean' in args
    if clean: args.remove('--clean')
    opts = {}
    for o in ('--roots', '--limit', '--rootlist'):
        if o in args:
            i = args.index(o); opts[o] = args[i + 1]; del args[i:i + 2]
    fam_name, table, outp, nproc = args[0], args[1], args[2], int(args[3]); gens = args[4:]
    if clean and (os.path.exists(outp) or opts):
        raise RuntimeError('--clean requires a new output and no partial-run options')
    if os.environ.get('MUB6_KERNEL'):
        raise RuntimeError('MUB6_KERNEL override forbidden in proof replay')
    log = lambda *a: print(time.strftime('%H:%M:%S'), *a, flush=True)
    log(f'checker: family {fam_name}; generator files {gens}')
    tsha = sha256(table); log(f'table {table} sha256 {tsha}')
    ck.load_table(table, verify=True, log=log)
    fam = ck.Fam(fam_name); fam.install()
    replay_context = rb.context(fam, table)
    recs, nlines = load_records(gens)
    log(f'{nlines} generator records read; {len(recs)} distinct certified boxes with hulls')
    need, missing, kc = walk(fam, recs)
    log(f'cover tree (checker selector): kept boxes per level {kc}; {len(need)} boxes to verify; {len(missing)} leaves without record')
    # resume: previously verified boxes in OUT
    done = rb.read_cache(outp, replay_context, recs)
    groups = {}
    for key in need:
        if key in done: continue
        lev, idx = key
        root = tuple(int(t) >> (lev - LROOT) for t in idx)
        groups.setdefault(root, []).append((lev, idx, recs[key][0], recs[key][1]))
    roots = sorted(groups)
    if '--roots' in opts:
        k, mm = opts['--roots'].split('/'); roots = [r for i, r in enumerate(roots) if i % int(mm) == int(k)]
    if '--rootlist' in opts:                       # restrict to the roots listed in a JSON file (work splitting)
        allowed = {tuple(r) for r in json.load(open(opts['--rootlist']))}
        roots = [r for r in roots if r in allowed]
    if '--limit' in opts: roots = roots[:int(opts['--limit'])]
    log(f'{len(done)} boxes already verified (resume); {len(roots)} roots / {sum(len(groups[r]) for r in roots)} boxes to do now')
    t0 = time.time(); nb = nok = 0
    ctx = get_context('fork') if sys.platform != 'darwin' else get_context('spawn')
    with ctx.Pool(nproc, initializer=init_worker, initargs=(fam_name, table, tsha, replay_context)) as pool, open(outp, 'a') as f:
        for i, res in enumerate(pool.imap_unordered(process_root, [(r, groups[r]) for r in roots])):
            f.write(json.dumps(res) + '\n'); f.flush()
            nb += len(res['boxes']); nok += sum(b['ok'] for b in res['boxes'])
            for b in res['boxes']:
                done[(b['lev'], tuple(b['idx']))] = b['ok']
                if not b['ok']: log('NOT VERIFIED', b)
            log(f"[{i+1}/{len(roots)}] root {res['root']}: {len(res['boxes'])} boxes, {sum(b['ok'] for b in res['boxes'])} ok, {res['time']}s; total {nok}/{nb}; elapsed {time.time()-t0:.0f}s")
    verified = {k: v for k, v in done.items()}
    if rb.context(fam, table) != replay_context:
        raise RuntimeError('proof environment changed during replay')
    full, bad = evaluate(fam, verified)
    code = {n: sha256(os.path.join(HERE, n)) for n in ('kernel.c', 'ck.py', 'run_check.py') if os.path.exists(os.path.join(HERE, n))}
    summ = {'family': fam_name, 'replay_context': replay_context, 'clean_start': clean, 'table_sha256': tsha, 'checker_code_sha256': code, 'platform': platform.platform(), 'generator_files': [(g, sha256(g)) for g in gens],
            'records': nlines, 'certified_records': len(recs), 'kept_per_level': kc, 'boxes_needed': len(need),
            'leaves_without_record': len(missing), 'verified_ok': int(sum(1 for k in need if done.get(k) is True)),
            'verified_fail': int(sum(1 for k in need if done.get(k) is False)), 'not_yet_checked': int(sum(1 for k in need if k not in done)),
            'cover_complete': bool(full), 'uncovered_level9_ancestors': len(bad), 'uncovered_examples': bad[:20]}
    log('SUMMARY', json.dumps(summ))
    with open(outp.replace('.jsonl', '_summary.json'), 'w') as f: json.dump(summ, f, indent=1)
    if not full: raise SystemExit(2)

if __name__ == '__main__':
    main()
