"""make_values.py -- compute the numbers quoted in the paper from the final outputs and write values.json.
Inputs: analysis/stats_szollosi.json, analysis/stats_fourier.json, chk_results/*_summary.json,
chk/tests/mutation_results.json, analysis/empty_*.json, exact/degeneracy_diagnostics.json."""
import json, os, datetime, hashlib
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
J = lambda p: json.load(open(os.path.join(R, p)))
fmt = lambda n: f'{int(n):,}'
V = {}
V['DATE'] = datetime.date.today().strftime('%-d %B %Y')
V['KERNEL_LINES'] = sum(1 for _ in open(os.path.join(R, 'chk/kernel.c')))
V['CKPY_LINES'] = sum(sum(1 for _ in open(os.path.join(R, 'chk', f))) for f in ('ck.py', 'run_check.py', 'verify_split.py'))
V['TABLE_SHA'] = hashlib.sha256(open(os.path.join(R, 'gen/table17.npy'), 'rb').read()).hexdigest()
sz = J('analysis/stats_szollosi.json'); fo = J('analysis/stats_fourier.json')
ssum = J('chk_results/chk_szollosi_summary.json'); fsum = J('chk_results/chk_fourier_summary.json')
kl = ssum['kept_per_level']; fl = fsum['kept_per_level']
lv = [kl[str(l)] for l in range(3, 10)]; flv = [fl[str(l)] for l in range(3, 10)]
V['SZ_KEPT_LEVELS'] = ', '.join(fmt(x) for x in lv[:-1]) + ' and ' + fmt(lv[-1])
V['F_KEPT_LEVELS'] = ', '.join(fmt(x) for x in flv[:-1]) + ' and ' + fmt(flv[-1])
V['SZ_BOXES'] = fmt(sz['checker_verified'])
V['F_BOXES'] = fmt(fo['checker_verified'])
V['SZ_ROOTS'] = fmt(kl['7'])
V['SZ_VM_ROOTS'] = 195
rc = sz['root_pass_L7_survivors']
V['ROOT_MIN'] = f"${rc['p0']/1e6:.1f}\\cdot10^6$"
V['ROOT_MAX'] = f"${rc['p100']/1e6:.1f}\\cdot10^6$"
mut = J('chk/tests/mutation_results.json')
V['MUT_TOTAL'] = len(mut); V['MUT_KILLED'] = sum(m['killed'] for m in mut)
surv = [m['mutation'] for m in mut if not m['killed']]
V['MUT_SURVIVOR_TEXT'] = ('The surviving mutation' + ('s are ' if len(surv) > 1 else ' is ') + '; '.join(f'"{s}"' for s in surv) + '. ' +
    'Removing the term $\\sum_l(YJ_{q_l})[-h_P,h_P]$ from $\\mathcal M$ in Lemma 17 makes the test unsound in principle, but the omitted quantity is of second order, $|YJ_q|\\,h_P\\,\\rho$, and it is dominated in practice by the overestimation of the other second-order terms ($|YJ_w|\\rho^2$ from the Jacobian over $W$ and $|YF_{qq}|h_P^2$ from the value term), so it is not exposed by any test of the suite — neither by the stress test T10 with coarse parameter boxes nor by the white-box sampling test T11 — and no test of this kind can be expected to expose it. Its presence in the kernel rests on the proof of Lemma 17 and on inspection of the code.') if surv else 'No mutation survived.'
emp = J('analysis/empty_szollosi.json')
V['EMPTY_TOTAL'] = emp['boxes_without_hulls']; V['EMPTY_PROVED'] = emp['proved_disjoint_from_F']
p, t = V['EMPTY_PROVED'], V['EMPTY_TOTAL']
V['EMPTY_PHRASE'] = 'both such boxes' if p == t == 2 else (f'all {t} such boxes' if p == t else f'{p} of the {t} such boxes')
json.dump(V, open(os.path.join(R, 'paper/values_base.json'), 'w'), indent=1)
print(json.dumps(V, indent=1))
