"""special_boxes.py -- checker results for the level-9 boxes containing alpha = 0, the exact ray point and alpha = 1.
Usage: python3 special_boxes.py OUT.json CHK1.jsonl [...] --gen G1.jsonl [...]"""
import sys, json
SPECIAL = {'alpha = 0': [0, 170, 0, 170], 'ray point 1/2 + i/(2 sqrt3)': [19, 137, 203, 322], 'alpha = 1 (Dita)': [0, 128, 0, 256]}
a = sys.argv[2:]; gi = a.index('--gen'); chks = a[:gi]; gens = a[gi + 1:]
boxes = {}
for p in chks:
    for line in open(p):
        r = json.loads(line)
        for b in r['boxes']: boxes[(b['lev'], tuple(b['idx']))] = b
gen = {}
for p in gens:
    for line in open(p):
        r = json.loads(line)
        if r.get('status') == 'certified': gen[(r['lev'], tuple(r['idx']))] = r
out = {}
for name, idx in SPECIAL.items():
    k = (9, tuple(idx)); b = boxes.get(k); g = gen.get(k)
    out[name] = {'idx': idx, 'checker': b, 'generator': ({kk: v for kk, v in g.items() if kk not in ('hull_lo', 'hull_hi')} if g else None)}
json.dump(out, open(sys.argv[1], 'w'), indent=1)
for n, v in out.items():
    b = v['checker']; g = v['generator']
    print(n, v['idx'], 'verified' if b and b['ok'] else 'NOT FOUND/NOT VERIFIED', b and b['nh'], b and b['graph_info'], b and b['cov_stats'], g and g.get('nsol'), g and g.get('max_hw'))
