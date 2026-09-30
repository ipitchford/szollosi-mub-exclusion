"""make_values2.py -- Section 9 numbers and tables (run after make_values.py); writes paper/values.json."""
import json, os, re, ast
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = lambda p: json.load(open(os.path.join(R, p)))
V = J('paper/values_base.json')
fmt = lambda n: f'{int(n):,}'

def log_cpu(paths):
    """sum of per-root times printed by the driver (each root is handled by one worker)"""
    tot = 0.0; n = 0
    for p in paths:
        if not os.path.exists(p): continue
        for line in open(p):
            m = re.match(r"\[\d+/\d+\] (\{.*\})\s+elapsed", line)
            if not m: continue
            try: d = ast.literal_eval(m.group(1))
            except Exception: continue
            tot += d.get('time', 0.0); n += 1
    return tot, n

sz = J('analysis/stats_szollosi.json'); fo = J('analysis/stats_fourier.json')
gl = lambda names: [os.path.join(R, 'logs', n) for n in names]
sz_gen_cpu, _ = log_cpu(gl(['gen_vm.log', 'gen_vm_tail.log', 'gen_vm_tail2.log', 'gen_mac.log', 'gen_mac2.log', 'gen_mac3.log']))
f_gen_cpu, _ = log_cpu(gl(['gen_fourier.log']))
sz_chk_cpu = sz['checker_time_roots_s']; f_chk_cpu = fo['checker_time_roots_s']
V['SZ_GEN_HOURS'] = f'{sz_gen_cpu / 3600:.0f}'
V['SZ_CHK_HOURS'] = f'{sz_chk_cpu / 3600:.0f}'
st = sz['generator_record_states']
V['SZ_FAILED'] = fmt(st.get('failed', 0))
dl = sz['checker_deepest_phase_level']
V['SZ_DEEPEST'] = f"{max(int(k) for k in dl)} (the deepest level needed was " + ', '.join(f"{k} for {fmt(v)} box" + ('' if v == 1 else 'es') for k, v in sorted(dl.items(), key=lambda t: int(t[0]))) + ')'
nsp = len(sz['checker_parameter_split'])
V['SZ_CHK_VERIFIED'] = (f"all {fmt(sz['checker_verified'])} certificates, {nsp} of them by parameter subdivision (Section 6.5) at parameter levels " +
                        f"{min(x[1] for x in sz['checker_parameter_split'])} to {max(x[1] for x in sz['checker_parameter_split'])}") if nsp else f"all {fmt(sz['checker_verified'])} certificates"
# parameter-split paragraph (from the split records and the numerical probe log)
import re as _re
splits = sz['checker_parameter_split']
probe = open(os.path.join(R, 'analysis/probe_split_boxes.log')).read() if os.path.exists(os.path.join(R, 'analysis/probe_split_boxes.log')) else ''
pb = _re.findall(r'== box ([0-9,]+)\n(?:.*\n)*?total (\d+) outside (\d+)', probe)
npar = len(_re.findall(r'^p=', probe, _re.M))
pb = [x for x in pb if int(x[1]) > 0]                  # boxes in which surface parameters were found
tot_s = sum(int(t) for _, t, _ in pb); tot_out = sum(int(o) for _, _, o in pb)
vac = [x for x in splits if x[2] == 0]
V['SZ_SPLIT_TEXT'] = (f"The {len(splits)} boxes that needed parameter subdivision all contain parameters with $z=w=1$ ($\\theta_x+\\theta_y\\equiv\\theta_u+\\theta_v\\equiv0$), where $\\alpha=1+2\\cos\\theta_x$ is real: they lie on the edge $\\arg\\alpha=0$ of the sector, with $\\alpha$ between about $0.16$ and the Diţă point $\\alpha=1$, and $\\theta_z,\\theta_w$ sit at the wrap-around point $0\\equiv2\\pi$ of the ordering conditions. "
                      + (f"For {len(vac)} of them no kept sub-box exists at the next level (the box does not meet $\\mathcal F$), so the verification is vacuous. " if vac else '')
                      + f"For the others the difficulty lies with the checker's first-order parametric Krawczyk test on the full parameter box, not with the hulls: a numerical probe, which is not part of the proof, sampled parameters on $\\Sigma$ in {len(pb)} of these boxes and found " + (f"all {fmt(tot_s)} numerically computed unbiased vectors inside the stored hulls" if tot_out == 0 else f"{tot_out} of {fmt(tot_s)} numerically computed unbiased vectors outside the stored hulls") + " (`analysis/probe_split_boxes.log`).")
V['SZ_ROOTS_WITH_CHILDREN'] = fmt(sz.get('roots_with_children', 0))
V['SZ_GEN_OUTCOME'] = (f"The generator certified every one of the {fmt(sz['certified_boxes'])} level-$9$ boxes at level $9$ "
                       f"({fmt(sz['second_attempt_Lmax10'])} of them in the second attempt with $L_{{\\max}}=10$); the only failed records are those of the incident described in Section 9.4.")
nh = sz['nsol_histogram']
hist = '; '.join(f"{k}: {fmt(v)}" for k, v in sorted(nh.items(), key=lambda t: int(t[0])))
hp = sz['hulls_per_box']; hw = sz['max_hull_halfwidth']; gc = sz['checker_graph']
rows = [
 ('certified level-9 boxes (generator), replayed and verified (checker)', f"{fmt(sz['certified_boxes'])} / {fmt(sz['checker_verified'])}"),
 ('numerically located solutions per box (complete histogram, value: boxes)', hist),
 ('boxes with rescued solutions; solutions rescued', f"{fmt(sz['boxes_with_rescue'])}; {fmt(sz['rescued_solutions'])}"),
 ('boxes with extra hulls from left-over clusters; extra hulls', f"{fmt(sz['boxes_with_extra_hulls'])}; {fmt(sz['extra_hulls_total'])}"),
 ('hulls per box: minimum, median, maximum', f"{hp['p0']:.0f}, {hp['p50']:.0f}, {hp['p100']:.0f}"),
 ('boxes without hulls', fmt(sz['boxes_without_hulls'])),
 ('largest hull half-width per box: median, 99th percentile, maximum', f"{hw['p50']:.3f}, {hw['p99']:.3f}, {hw['p100']:.3f}"),
 ('six-cliques of $O$ (checker): median and maximum per box; total', f"{gc['cliques']['p50']:.0f}, {gc['cliques']['p100']:.0f}; {fmt(gc['cliques_total'])}"),
 ('edges of $O$ and $B$ (checker), maxima', f"{fmt(gc['edgesO_max'])}, {fmt(gc['edgesB_max'])}"),
 ('root-pass survivors at level 7 (checker): minimum, median, maximum', f"{fmt(sz['root_pass_L7_survivors']['p0'])}, {fmt(sz['root_pass_L7_survivors']['p50'])}, {fmt(sz['root_pass_L7_survivors']['p100'])}"),
 ('CPU time: generator, checker (hours)', f"{sz_gen_cpu / 3600:.1f}, {sz_chk_cpu / 3600:.1f}"),
]
V['SZ_TABLE'] = '\n'.join(f'| {a} | {b} |' for a, b in rows)
fh = fo['nsol_histogram']; fhp = fo['hulls_per_box']; fhw = fo['max_hull_halfwidth']; fgc = fo['checker_graph']
frows = [
 ('certified level-9 boxes (generator), replayed and verified (checker)', f"{fmt(fo['certified_boxes'])} / {fmt(fo['checker_verified'])}"),
 ('numerically located solutions per box', '; '.join(f"{k}: {fmt(v)}" for k, v in sorted(fh.items(), key=lambda t: int(t[0])))),
 ('hulls per box: minimum, median, maximum', f"{fhp['p0']:.0f}, {fhp['p50']:.0f}, {fhp['p100']:.0f}"),
 ('largest hull half-width per box: median, maximum', f"{fhw['p50']:.3f}, {fhw['p100']:.3f}"),
 ('six-cliques of $O$ (checker): median and maximum per box', f"{fgc['cliques']['p50']:.0f}, {fgc['cliques']['p100']:.0f}"),
 ('deepest phase level in the checker (level: boxes)', '; '.join(f"{k}: {fmt(v)}" for k, v in sorted(fo['checker_deepest_phase_level'].items(), key=lambda t: int(t[0])))),
 ('CPU time: generator, checker (hours)', f"{f_gen_cpu / 3600:.2f}, {f_chk_cpu / 3600:.2f}"),
]
V['F_TABLE'] = '\n'.join(f'| {a} | {b} |' for a, b in frows)
fsum = J('chk_results/chk_fourier_summary.json')
cs = {}
for line in open(os.path.join(R, 'chk_results/chk_fourier.jsonl')):
    r = json.loads(line)
    for bx in r['boxes']: cs[tuple(bx['idx'])] = bx['graph_info'].get('cliques', 0)
n8 = sum(1 for v in cs.values() if v == 8)
c00 = cs.get((0, 0)); c16 = cs.get((85, 0)); cmax = max(cs.items(), key=lambda t: t[1])
V['F_TEXT'] = (f"The selector of Section 4.3 keeps {V['F_KEPT_LEVELS']} boxes at levels $3,\\dots,9$. The generator certified all "
               f"{fmt(fo['certified_boxes'])} level-$9$ boxes at level $9$, and the checker verified all of them; its cover evaluation is complete "
               f"(`chk_results/chk_fourier_summary.json`). By Theorem 7 and Proposition 12 this proves Theorem 3.\n\n"
               f"The certificates agree with the numerical observations of Jaming et al., which they do not depend on. In every box the generator located "
               f"{min(fh, key=int)} unbiased vectors numerically (their Numerical Evidence 1.3: $48$ for all $(a,b)$). The six-cliques of the graph $O$ are the candidate orthonormal bases: "
               f"in {fmt(n8)} of the {fmt(len(cs))} boxes there are exactly $8$ of them (their observation for generic parameters), in the box at $(a,b)=(0,0)$ there are {c00} "
               f"(they report $16$ bases for $F_6$), and in the box containing $(a,b)=(1/6,0)$ there are exactly {c16} (they report $70$). The largest number, {cmax[1]}, occurs in the "
               f"neighbouring box {cmax[0]}, where hulls of nearby vectors overlap.")
sp = J('analysis/special_szollosi.json')
parts = []
NAMES = {'alpha = 0': '$\\alpha=0$', 'ray point 1/2 + i/(2 sqrt3)': 'the ray point $\\alpha=\\tfrac12+\\tfrac{i}{2\\sqrt3}$', 'alpha = 1 (Dita)': 'the Diţă point $\\alpha=1$'}
for name, v in sp.items():
    b = v['checker']; g = v['generator']
    if not b: parts.append(f"{NAMES.get(name, name)}: box {v['idx']} not found"); continue
    how = "verified directly" if not b.get('split_detail') else f"verified by parameter subdivision to level {b['split_detail']['split_level']}"
    parts.append(f"the box {tuple(v['idx'])} containing {NAMES.get(name, name)} was {how}, with {b['nh']} hulls ({g.get('nsol') if g else '?'} numerical solutions, largest hull half-width {g.get('max_hw', 0):.3f}); "
                 f"its graph $O$ has {b['graph_info'].get('edgesO')} edges and {b['graph_info'].get('cliques')} six-clique{'' if b['graph_info'].get('cliques') == 1 else 's'} and $B$ has {fmt(b['graph_info'].get('edgesB'))} edges")
V['SZ_SPECIAL'] = 'At the special parameters: ' + '; '.join(parts) + '. In each case no six-clique of $O$ has a partner clique across $B$.'
V['SZ_SPLIT_COUNT'] = len(sz['checker_parameter_split'])
json.dump(V, open(os.path.join(R, 'paper/values.json'), 'w'), indent=1)
print(json.dumps({k: (v if len(str(v)) < 300 else str(v)[:300] + '...') for k, v in V.items()}, indent=1))
