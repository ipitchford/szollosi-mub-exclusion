"""
verify_split.py -- CHECKER: verification of a certificate by parameter subdivision.

A certificate (hull family H_P) for a parameter box P of level l claims that U(H(p)) is covered by H_P for every
p in the half-open ownership cell Q intersected with the fundamental domain F. By Lemma 20,
every such p has a retained descendant owner within Q at every finer level. Hence it suffices to verify coverage
of the SAME hull family on every kept sub-box of P at level l + d (checker selector), which is often possible when
the verification on P itself fails because the parameter box is too wide for the exclusion tests near a
degenerate parameter.  The combinatorial conditions (Sep) and (Cl) concern the hulls only and are checked once.

Usage: python3 verify_split.py FAMILY TABLE.npy OUT.jsonl GEN1.jsonl [...] --boxes BOXES.json [--depth 2]
BOXES.json: list of [lev, idx].  Output records have the same format as run_check.py (one root per line), with
'split' information; a later ok=True record for a box supersedes an earlier failed one in the cover evaluation.
"""
import sys, os, json, time, hashlib, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ck
import receipt_binding as rb
from run_check import load_records, sha256, LROOT

def main():
    a = sys.argv[1:]
    depth = 2
    if '--depth' in a:
        i = a.index('--depth'); depth = int(a[i + 1]); del a[i:i + 2]
    i = a.index('--boxes'); boxes = json.load(open(a[i + 1])); del a[i:i + 2]
    fam_name, table, outp = a[0], a[1], a[2]; gens = a[3:]
    if os.environ.get('MUB6_KERNEL'):
        raise RuntimeError('MUB6_KERNEL override forbidden in proof replay')
    if depth < 1 or any(lev + depth > 12 for lev, idx in boxes):
        raise ValueError('subdivision depth outside certified levels')
    log = lambda *x: print(time.strftime('%H:%M:%S'), *x, flush=True)
    tsha = sha256(table); log('table sha256', tsha)
    ck.load_table(table, verify=True, log=log)
    fam = ck.Fam(fam_name); fam.install(); sel = ck.selector(fam); m = fam.m
    recs, n = load_records(gens)
    replay_context = rb.context(fam, table)
    rb.read_cache(outp, replay_context, recs)
    offs = np.indices((2,) * m).reshape(m, -1).T
    for lev, idx in boxes:
        key = (lev, tuple(idx)); hlo, hhi = recs[key]
        hlo = np.asarray(hlo, dtype=np.int64).reshape(-1, 5); hhi = np.asarray(hhi, dtype=np.int64).reshape(-1, 5)
        hsha = hashlib.sha256(json.dumps([hlo.tolist(), hhi.tolist()]).encode()).hexdigest()[:16]
        root = tuple(int(t) >> (lev - LROOT) for t in idx)
        t0 = time.time()
        pc, ph = ck.pbox(np.array(root), LROOT)
        levels, counts = ck.root_pass(pc, ph)
        ok_g, gi = ck.verify_graphs(hlo, hhi)
        result = None
        for d in range(1, depth + 1):
            cur = np.asarray(idx, dtype=np.int64)[None, :]
            for k in range(d):
                ch = (cur[:, None, :] * 2 + offs[None]).reshape(-1, m)
                cur = ch[sel(ch, lev + k + 1)]
            stats = []; all_ok = True
            for c in cur:
                okc, st = ck.verify_coverage(fam, c, lev + d, hlo, hhi, levels, L7=LROOT)
                stats.append({'idx': c.tolist(), 'ok': okc, 'cov_stats': st}); all_ok &= okc
                if not okc:
                    log(f'  sub-box {c.tolist()} (level {lev + d}) not verified: {st}')
                    break
            log(f'box {lev} {idx}: split to level {lev + d}, {len(cur)} kept sub-boxes, coverage {all_ok}')
            if all_ok:
                result = {'split_level': lev + d, 'sub_boxes': len(cur), 'sub_stats': stats}; break
        rec = {'lev': lev, 'idx': list(idx), 'root': list(root), 'hulls_sha': hsha, 'ok': bool(result is not None and ok_g),
               'coverage': result is not None, 'graphs': ok_g, 'nh': int(len(hlo)),
               'cov_stats': {'split_level': result['split_level'] if result else None, 'sub_boxes': result['sub_boxes'] if result else None},
               'split_detail': result, 'graph_info': gi, 'time': round(time.time() - t0, 1),
               'input_sha256': rb.binding(replay_context, key, [hlo.tolist(), hhi.tolist()])}
        if rb.context(fam, table) != replay_context:
            raise RuntimeError('proof environment changed during subdivision')
        with open(outp, 'a') as f:
            f.write(json.dumps({'root': list(root), 'root_counts': counts, 'boxes': [rec], 'time': rec['time'], 'mode': 'parameter-split', 'replay_context': replay_context}, default=int) + '\n')
        log('RESULT', rec['lev'], rec['idx'], 'ok' if rec['ok'] else 'NOT VERIFIED')

if __name__ == '__main__':
    main()
