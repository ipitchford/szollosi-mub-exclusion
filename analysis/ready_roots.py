"""ready_roots.py -- work splitting helper (not part of the proof): list the level-7 roots (checker selector) all of
whose kept level-9 children have a certified record in the given generator files.
Usage: python3 ready_roots.py FAMILY OUT.json GEN1.jsonl [GEN2.jsonl ...]"""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'chk'))
import ck
fam_name, out = sys.argv[1], sys.argv[2]
ck.load_table(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gen', 'table17.npy'), verify=False)
fam = ck.Fam(fam_name); fam.install()
have = set()
for p in sys.argv[3:]:
    for line in open(p):
        try: r = json.loads(line)
        except Exception: continue
        if r.get('status') == 'certified' and r.get('ok') is True: have.add((r['lev'], tuple(r['idx'])))
roots = ck.selected_boxes(fam, 7)
m = fam.m; offs = np.indices((4,) * m).reshape(m, -1).T
ready = []
for r in roots:
    ch = r * 4 + offs; ch = ch[ck.selector(fam)(ch, 9)]
    if len(ch) and all((9, tuple(int(t) for t in c)) in have for c in ch): ready.append([int(t) for t in r])
json.dump(ready, open(out, 'w'))
print(len(roots), 'roots,', len(ready), 'ready')
