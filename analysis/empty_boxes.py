"""
empty_boxes.py -- diagnostic for certified boxes with NO hulls (referee minor comment on empty hulls).
Such a box is certified because the checker's coverage verification excluded every phase box, i.e. no unbiased
vector exists for any p in the box on the constraint surface; no statement about P and Sigma is needed.
Separately, this script tries to prove P cap F = empty with the checker's selector by subdivision (to level 14).
Usage: python3 empty_boxes.py FAMILY OUT.json CHK1.jsonl [CHK2.jsonl ...]
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'chk'))
import ck, numpy as np

def main():
    fam_name, out = sys.argv[1], sys.argv[2]
    ck.load_table(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'gen', 'table17.npy'), verify=False)
    fam = ck.Fam(fam_name); fam.install()
    last = {}
    for p in sys.argv[3:]:
        for line in open(p):
            r = json.loads(line)
            for b in r['boxes']: last[(b['lev'], tuple(b['idx']))] = b
    empty = [(k[0], list(k[1])) for k, b in last.items() if b['ok'] and b['nh'] == 0]
    res = []
    for lev, idx in empty:
        ok, lv = ck.prove_empty(fam, np.array(idx), lev, max_lev=14)
        res.append({'lev': lev, 'idx': idx, 'proved_disjoint_from_F': bool(ok), 'level_reached': int(lv)})
    S = {'family': fam_name, 'boxes_without_hulls': len(empty), 'proved_disjoint_from_F': sum(r['proved_disjoint_from_F'] for r in res), 'boxes': res}
    json.dump(S, open(out, 'w'), indent=1)
    print(json.dumps({k: v for k, v in S.items() if k != 'boxes'}))

if __name__ == '__main__':
    main()
