"""V3 regression tests: binding, canonical owners, and exact mixed-term fixture."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from fractions import Fraction as F
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import ck
import receipt_binding as rb
import numpy as np

def require(ok, message):
    if not ok: raise RuntimeError(message)

def rejected(fn):
    try: fn()
    except (RuntimeError, ValueError): return
    raise RuntimeError('negative control accepted')

root = Path(__file__).resolve().parents[2]
ck.load_table(root / 'gen/table17.npy', verify=True)
fam = ck.Fam('fourier'); fam.install()
ctx = rb.context(fam, root / 'gen/table17.npy')
key = (9, (0, 0)); hulls = ([[0]*5], [[1]*5]); recs = {key: hulls}
with tempfile.TemporaryDirectory() as td:
    path = Path(td)/'cache.jsonl'
    box = {'lev': 9, 'idx': [0,0], 'ok': True, 'input_sha256': rb.binding(ctx,key,hulls)}
    row = {'replay_context':ctx, 'boxes':[box]}
    path.write_text(json.dumps(row)+'\n')
    require(rb.read_cache(path,ctx,recs)[key], 'same-input resume failed')
    rejected(lambda:rb.read_cache(path,ctx,{key:([],[])}))
    for field in ('family','table_sha256','kernel_sha256','code_sha256','family_data','numpy'):
        changed=copy.deepcopy(ctx); changed[field]='changed'
        rejected(lambda:rb.read_cache(path,changed,recs))
    path.write_text(json.dumps({'boxes':[box]})+'\n')
    rejected(lambda:rb.read_cache(path,ctx,recs))
    # Exercise the real driver with the referee's changed Fourier origin and legacy cache.
    cert=Path(td)/'altered.jsonl'; cache=Path(td)/'legacy.jsonl'
    with cert.open('w') as out:
        for line in (root/'results/gen_fourier.jsonl').open():
            rec=json.loads(line)
            if rec['lev']==9 and rec['idx']==[0,0]: rec['hull_lo']=[];rec['hull_hi']=[]
            out.write(json.dumps(rec)+'\n')
    cache.write_bytes((root/'chk_results/chk_fourier.jsonl').read_bytes())
    cmd=[sys.executable,str(root/'chk/run_check.py'),'fourier',str(root/'gen/table17.npy'),str(cache),'1',str(cert)]
    p=subprocess.run(cmd,capture_output=True,text=True)
    require(p.returncode!=0 and 'unbound or mismatched replay context' in p.stderr,'driver accepted stale cache')
    require(not cache.with_name('legacy_summary.json').exists(),'failure emitted a success summary')
    p=subprocess.run(cmd+['--clean'],capture_output=True,text=True)
    require(p.returncode!=0 and '--clean requires' in p.stderr,'clean mode accepted existing output')
    p=subprocess.run(cmd,env=os.environ|{'MUB6_KERNEL':'untrusted'},capture_output=True,text=True)
    require(p.returncode!=0 and 'override forbidden' in p.stderr,'proof accepted kernel override')

fam=ck.Fam('szollosi');fam.install()
require(ck.select_szollosi(np.array([[34,106,264,263],[34,106,264,264]]),9).tolist()==[False,True], 'boundary ownership regression')
# Exact generic Krawczyk assembly; does not claim to kill a kernel-specific mutant.
w,q=F(5,6),F(2); require((1+q)*w-F(5,2)==0,'fixture root')
bad=F(5,4); lo,hi=bad-F(1,2),bad+F(1,2)
require(bad>1 and lo<=w<=hi and -1<=w<=1,'mixed derivative is necessary')
print('PASS: content binding, legacy-cache driver rejection, clean mode, override rejection, boundary ownership, mixed-term assembly')
