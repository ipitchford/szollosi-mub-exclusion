"""make_figure.py -- draw Figure 1 (data flow and trust) with matplotlib."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

W, H = 8.6, 5.5
fig = plt.figure(figsize=(6.4, 6.4 * H / W))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(0, H); ax.axis('off')
FS = 7.2

def cluster(x0, y0, x1, y1, fc, label):
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc=fc, ec='#888888', lw=0.8))
    ax.text(x0 + 0.08, y1 - 0.08, label, ha='left', va='top', fontsize=FS, style='italic', color='#333333')

def node(x, y, w, h, text, fc='white', bold_first=True):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.02,rounding_size=0.06',
                                fc=fc, ec='#555555', lw=0.8))
    lines = text.split('\n')
    ax.text(x, y, '\n'.join(lines), ha='center', va='center', fontsize=FS, linespacing=1.25,
            family='DejaVu Sans')
    return (x, y, w, h)

def arrow(p, q, label=None, lpos=None, ha='center'):
    ax.annotate('', xy=q, xytext=p, arrowprops=dict(arrowstyle='-|>', color='#333333', lw=0.8,
                                                    shrinkA=0, shrinkB=0, mutation_scale=8))
    if label:
        ax.text(*lpos, label, fontsize=FS - 0.4, ha=ha, va='center', style='italic', color='#222222')

cluster(1.30, 4.00, 6.05, 5.42, '#fbf3d0', 'Untrusted: may be wrong without affecting the theorems')
cluster(1.30, 0.08, 6.05, 3.70, '#e2ebf8', 'Trusted: the proof depends on these')

node(0.62, 3.85, 1.05, 0.62, 'gen/\ntable17.npy\nunit table', fc='#f2f2f2')
node(2.20, 4.55, 1.60, 0.72, 'gen/core.py\nradius heuristics,\nArb table check')
node(3.98, 4.55, 1.74, 0.72, 'gen/ generator:\nNewton iteration,\nbranch and bound, hulls')
node(5.47, 4.55, 1.00, 0.72, 'results/\ncertificates:\nboxes, hulls')

node(2.20, 2.95, 1.60, 0.62, 'mpmath intervals:\nunit table and\nconstants')
node(2.20, 1.85, 1.60, 1.05, 'chk/kernel.c\none-ulp outward\nkernel: modulus,\nmean-value, Krawczyk,\ncontainment, pairs')
node(4.70, 2.35, 2.10, 1.25, 'chk/ck.py, run_check.py,\nverify_split.py:\nfamilies from formulas,\nselector, coverage branch\nand bound, cliques, cover tree')
node(2.20, 0.62, 1.60, 0.82, 'exact/, fourier/\nSymPy and exact\nrationals (Lemma 9,\nProps 10, 12)')
node(4.70, 0.62, 2.10, 0.82, 'Mathematics of Sections 2–6\n(Lemmas 4–20; Props 6, 10,\n12, 18; Theorem 7)')

cx, cy, hw, hh = 7.35, 2.05, 1.05, 0.78
ax.add_patch(Polygon([(cx - hw, cy), (cx, cy + hh), (cx + hw, cy), (cx, cy - hh)], closed=True,
                     fc='#eeeeee', ec='#555555', lw=0.8))
ax.text(cx, cy, 'every certificate\nverified and\ncover complete', ha='center', va='center', fontsize=FS, linespacing=1.25)
node(7.35, 0.55, 1.90, 0.62, 'Theorem 1 (Szöllősi)\nTheorem 3 (Fourier)', fc='#f2f2f2')

arrow((0.62, 4.16), (1.40, 4.50))
arrow((0.62, 3.54), (1.40, 2.98), 'verified\nindependently', (0.10, 3.05), ha='left')
arrow((3.00, 4.55), (3.11, 4.55))
arrow((4.85, 4.55), (4.97, 4.55))
arrow((5.47, 4.19), (5.30, 2.975), 'read as data only', (5.52, 3.85), ha='left')
arrow((3.00, 2.95), (3.65, 2.75))
arrow((3.00, 1.85), (3.65, 2.05))
arrow((3.00, 0.62), (3.65, 0.62))
arrow((5.75, 2.35), (cx - hw + 0.02, cy + 0.02))
arrow((5.75, 0.80), (cx - 0.52, cy - hh / 2 - 0.02))
arrow((cx, cy - hh), (cx, 0.86))
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_dependency.png'), dpi=300)
fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_dependency.pdf'))
