"""Permanent interval-supported owner fixture adapted from supplied review.

No claim of a missing point: the discarded closed cell shares its upper face
with the retained canonical owner. Explicit checks also execute under -O.
"""
from pathlib import Path
import sys
import numpy as np
from mpmath import iv

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'chk'))
import ck

def require(value, label):
    if not value:
        raise RuntimeError(label)

ck.load_table(str(ROOT / 'gen/table17.npy'), verify=True)
family = ck.Fam('szollosi')
family.install()
iv.dps = 50
t = iv.pi / 32
a = 2 * iv.cos(t) - iv.cos(2*t)
b = 2 * iv.sin(t) + iv.sin(2*t)
require(a.a > 0 and b.a > 0 and (a-iv.sqrt(3)*b).a > 0, 'sector')

def g(k):
    theta = 2 * iv.pi * k / 512
    return iv.sin(3*theta/2) - a*iv.sin(theta/2) - b*iv.cos(theta/2)

for lo in [34, 106]:
    left, right = g(lo), g(lo+1)
    require((left.b < 0 and right.a > 0) or
            (left.a > 0 and right.b < 0), 'root bracket')
# IVT gives two distinct unit roots; Vieta and xyz=1 give the third.
require(1-142/512 > 107/512, 'third root ordering')
indices = np.array([[34,106,264,263], [34,106,264,264]], dtype=np.int64)
require(ck.select_szollosi(indices, 9).tolist() == [False, True], 'owner')
print('PASS: interval sector, two root brackets, third-root ordering, owner retention')
