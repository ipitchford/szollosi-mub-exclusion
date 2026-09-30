"""
mutate.py -- mutation testing of the checker kernel: each mutation removes or weakens one soundness-relevant
ingredient; the mutated kernel is compiled separately and the test-suite (test_kernel.py) must FAIL on it.
A mutation that survives would show a test gap.  The trusted kernel.c itself is never modified.
"""
import os, sys, subprocess, tempfile, json, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); SRC = os.path.join(HERE, '..', 'kernel.c')
MUTATIONS = [
 ('no outward rounding (up1 = identity)', 'static inline double up1(double x) {\n    if (x != x) return x;',
  'static inline double up1(double x) {\n    return x;\n    if (x != x) return x;'),
 ('table error ignored', 'TR[m].re.lo = dn1(re - TR_EPS); TR[m].re.hi = up1(re + TR_EPS);', 'TR[m].re.lo = re; TR[m].re.hi = re;'),
 ('sine bound used beyond its range', 'double sm = (w <= NG / 4) ? TR[w].im.hi : 1.0;', 'double sm = TR[w & NGMASK].im.hi;'),
 ('modulus test without chord radius', 'double R = angle_up(W);', 'double R = 0.0 * angle_up(W);'),
 ('modulus test ignores parameter width', 'W += wp[k][j] + ((k == 0) ? 0 : n[k - 1]);', 'W += ((k == 0) ? 0 : n[k - 1]);'),
 ('mean-value test drops parameter derivatives', 'F = iv_add(F, iv_sym_mul(iv_scal(2.0, c_re_conj_mul(SP, d)), hq));', ';'),
 ('Krawczyk drops the parameter term of F', 'V = iv_add(V, iv_sym_mul(yf, hq));', ';'),
 ('Krawczyk drops the parameter term of J', 'E = iv_add(E, iv_sym_mul(yjq, hq));', ';'),
 ('Krawczyk uses the box radius without outward step', 'for (int k = 0; k < 5; k++) rad[k] = angle_up(nn[k]);', 'for (int k = 0; k < 5; k++) rad[k] = 0.5 * angle_up(nn[k]);'),
 ('containment off by one', 'if (off + (hi[5 * b + k] - lo[5 * b + k]) > L) ok = 0;', 'if (off + (hi[5 * b + k] - lo[5 * b + k]) > L + 1) ok = 0;'),
 ('pair bounds: orthogonality flag inverted', 'uint8_t o = (uint8_t)(p.lo <= 0.0)', 'uint8_t o = (uint8_t)(p.lo > 0.0)'),
 ('self-check too permissive', 'selfok[i] = (uint8_t)(S.re.lo > 0.0 && q.lo > 6.0);', 'selfok[i] = 1;'),
 ('selector without Taylor remainder', 'if (num > up1(nn * rem)) ec = 1;', 'if (num > 0.0) ec = 1;'),
 ('selector sector test with wrong radius', 'double r = up1(4.0 * hw);', 'double r = up1(1.0 * hw);'),
 ('pair bounds: unbiasedness flag dropped', 'bb = (uint8_t)(p.lo <= 6.0 && 6.0 <= p.hi);', 'bb = (uint8_t)(p.lo <= 6.0 && 6.0 <= p.hi && 0);'),
 ('mean-value test drops phase derivatives', 'F = iv_add(F, iv_sym_mul(iv_scal(2.0, c_re_conj_mul(SP, c_muli(tP[k]))), rad[k - 1]));', ';'),
 ('Krawczyk ignores the dependent-parameter width', 'for (int i = 0; i < nd; i++) rad[5 + i] = hq;', 'for (int i = 0; i < nd; i++) rad[5 + i] = 0.0;'),
]
def one(args):
    i, name, a, b, src, tmp = args
    m = src.replace(a, b, 1)
    cfile = os.path.join(tmp, 'k%d.c' % i); so = os.path.join(tmp, 'libmut_%d%s' % (i, '.dylib' if sys.platform == 'darwin' else '.so'))
    open(cfile, 'w').write(m)
    flags = ['-dynamiclib'] if sys.platform == 'darwin' else ['-shared']
    subprocess.run(['cc', '-O2', '-fPIC'] + flags + ['-ffp-contract=off', '-o', so, cfile, '-lm'], check=True, capture_output=True)
    env = dict(os.environ, MUB6_KERNEL=so)
    r = subprocess.run([sys.executable, os.path.join(HERE, 'test_kernel.py')], env=env, capture_output=True, text=True)
    failed = [l for l in r.stdout.splitlines() if l.startswith('FAIL')]
    tests = sorted({l.split()[1] for l in failed})
    return {'mutation': name, 'killed': r.returncode != 0, 'failing_tests': tests, 'n_failures': len(failed)}

def run(nproc=1):
    tmp = tempfile.mkdtemp()
    src = open(SRC).read()
    for name, a, b in MUTATIONS: assert a in src, name
    jobs = [(i, name, a, b, src, tmp) for i, (name, a, b) in enumerate(MUTATIONS)]
    if nproc > 1:
        from multiprocessing import Pool
        with Pool(nproc) as pool: out = pool.map(one, jobs)
    else:
        out = [one(j) for j in jobs]
    for o in out:
        print(f"{'KILLED ' if o['killed'] else 'SURVIVED'}  {o['mutation']}: failing tests {o['failing_tests']} ({o['n_failures']} failures)", flush=True)
    shutil.rmtree(tmp)
    json.dump(out, open(os.path.join(HERE, 'mutation_results.json'), 'w'), indent=1)
    print(f"{sum(o['killed'] for o in out)} of {len(out)} mutations killed")
if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
