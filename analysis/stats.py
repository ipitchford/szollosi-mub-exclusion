"""
stats.py -- statistics of the generator runs and the checker replay (for Section 9 of the paper).
Usage: python3 stats.py FAMILY OUT.json --gen G1.jsonl [G2.jsonl ...] --chk C1.jsonl [C2.jsonl ...]
Generator records: the record used is the LAST certified record of each box in the order of the files given
(the same rule as the checker's loader); records are matched with the checker's output by the hull digest.
"""
import sys, json, hashlib, numpy as np
from collections import Counter

def load_gen(paths):
    recs = {}; allrecs = []; status = Counter()
    for p in paths:
        for line in open(p):
            try: r = json.loads(line)
            except Exception: continue
            allrecs.append(r); status[r.get('status', 'none')] += 1
            if r.get('status') == 'certified' and r.get('ok') is True and 'hull_lo' in r:
                recs.setdefault((r['lev'], tuple(r['idx'])), []).append(r)
    return recs, allrecs, status

def load_chk(paths):
    boxes = {}; roots = []
    for p in paths:
        for line in open(p):
            try: r = json.loads(line)
            except Exception: continue
            roots.append(r)
            for b in r.get('boxes', []):
                boxes[(b['lev'], tuple(b['idx']))] = b
    return boxes, roots

def hsha(r):
    hlo = np.asarray(r['hull_lo'], dtype=np.int64).reshape(-1, 5); hhi = np.asarray(r['hull_hi'], dtype=np.int64).reshape(-1, 5)
    return hashlib.sha256(json.dumps([hlo.tolist(), hhi.tolist()]).encode()).hexdigest()[:16]

def q(a, ps=(0, 50, 90, 99, 100)):
    a = np.asarray(a, dtype=float)
    return {f'p{p}': float(np.percentile(a, p)) for p in ps} if len(a) else {}

def main():
    fam = sys.argv[1]; out = sys.argv[2]
    a = sys.argv[3:]
    gi = a.index('--gen'); ci = a.index('--chk')
    gens = a[gi + 1:ci]; chks = a[ci + 1:]
    recs, allrecs, status = load_gen(gens)
    boxes, roots = load_chk(chks)
    all_certified = len(recs)
    S_dups = sum(1 for v in recs.values() if len(v) > 1)
    chosen = {}; mism = []
    for k, b in boxes.items():                                  # the certificate actually replayed by the checker
        if k not in recs: continue
        cands = [r for r in recs[k] if hsha(r) == b.get('hulls_sha')] if b.get('hulls_sha') else recs[k][-1:]
        if cands: chosen[k] = cands[-1]
        else: mism.append(list(k)); chosen[k] = recs[k][-1]
    recs = chosen
    S = {'family': fam, 'generator_files': gens, 'checker_files': chks, 'certified_records_all_levels': all_certified}
    S['generator_record_states'] = dict(status)
    S['certified_boxes'] = len(recs)
    S['levels'] = dict(Counter(k[0] for k in recs))
    S['second_attempt_Lmax10'] = sum(1 for r in recs.values() if r.get('Lmax') == 10)
    S['nsol_histogram'] = dict(sorted(Counter(r.get('nsol') for r in recs.values()).items()))
    S['hulls_per_box'] = q([r['hulls'] for r in recs.values()])
    S['boxes_without_hulls'] = sum(1 for r in recs.values() if r['hulls'] == 0)
    S['max_hull_halfwidth'] = q([r.get('max_hw', 0) for r in recs.values() if r['hulls'] > 0])
    S['rescued_solutions'] = int(sum(r.get('rescued', 0) for r in recs.values()))
    S['boxes_with_rescue'] = sum(1 for r in recs.values() if r.get('rescued', 0) > 0)
    S['boxes_with_extra_hulls'] = sum(1 for r in recs.values() if r.get('extra_hulls', 0) > 0)
    S['extra_hulls_total'] = int(sum(r.get('extra_hulls', 0) for r in recs.values()))
    S['gen_cliques'] = q([r['cert'].get('cliques', 0) for r in recs.values() if 'cert' in r])
    S['gen_cliques_total'] = int(sum(r['cert'].get('cliques', 0) for r in recs.values() if 'cert' in r))
    S['gen_time_s_total'] = float(sum(r.get('time', 0) for r in allrecs))
    S['gen_time_per_box'] = q([r.get('time', 0) for r in recs.values()])
    S['survivors_L7_after_retest'] = q([r['L7'][1] for r in recs.values() if 'L7' in r])
    # checker
    verified = {k: b for k, b in boxes.items() if b['ok']}
    S['checker_boxes'] = len(boxes); S['checker_verified'] = len(verified)
    S['checker_failed'] = [list(k) for k, b in boxes.items() if not b['ok']][:50]
    S['boxes_with_several_generator_records'] = S_dups
    S['hull_digest_mismatches'] = mism[:20]; S['hull_digest_mismatch_count'] = len(mism)
    deepest = Counter(); maxlev = []
    for b in boxes.values():
        levs = [int(k[1:]) for k in b['cov_stats'] if k.startswith('L') and k[1:].isdigit()]
        if not levs and b.get('split_detail'):
            levs = [int(k[1:]) for s_ in b['split_detail']['sub_stats'] for k in s_['cov_stats'] if k.startswith('L') and k[1:].isdigit()]
        if levs: deepest[max(levs)] += 1
    S['checker_deepest_phase_level'] = dict(sorted(deepest.items()))
    S['checker_parameter_split'] = [[list(k), b['split_detail']['split_level'], b['split_detail']['sub_boxes']] for k, b in boxes.items() if b.get('split_detail')]
    S['checker_time_boxes_s'] = float(sum(b['time'] for b in boxes.values()))
    S['checker_time_roots_s'] = float(sum(r['time'] for r in roots))
    S['checker_roots'] = len(roots)
    S['roots_with_children'] = len({tuple(b['root']) for b in boxes.values()})
    rc = [r['root_counts'][-1][2] for r in roots if r.get('root_counts')]
    S['root_pass_L7_survivors'] = q(rc)
    S['checker_graph'] = {'cliques': q([b['graph_info'].get('cliques', 0) for b in boxes.values() if b.get('graph_info')]),
                          'cliques_total': int(sum(b['graph_info'].get('cliques', 0) for b in boxes.values() if b.get('graph_info'))),
                          'edgesO_max': int(max([b['graph_info'].get('edgesO', 0) for b in boxes.values() if b.get('graph_info')] or [0])),
                          'edgesB_max': int(max([b['graph_info'].get('edgesB', 0) for b in boxes.values() if b.get('graph_info')] or [0]))}
    json.dump(S, open(out, 'w'), indent=1)
    print(json.dumps(S, indent=1))

if __name__ == '__main__':
    main()
