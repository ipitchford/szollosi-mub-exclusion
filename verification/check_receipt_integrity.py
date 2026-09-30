"""Check recorded replay identity, not re-execute its mathematical computation."""
import hashlib, json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok,message):
    if not ok: raise RuntimeError(message)

for family in ('fourier','szollosi'):
    p=ROOT/'replay'/f'{family}_summary.json'
    s=json.loads(p.read_text());c=s['replay_context']
    require(s['family']==family and c['family']==family,'family identity')
    require(s['cover_complete'] is True,'cover not complete')
    require(s['verified_fail']==0 and s['not_yet_checked']==0,'unverified certificates')
    require(s['leaves_without_record']==0,'missing leaves')
    require(s['verified_ok']==s['boxes_needed'],'count mismatch')
    for name,h in c['code_sha256'].items():
        require(digest(ROOT/'chk'/name)==h,'changed checker '+name)
    require(digest(ROOT/'gen/table17.npy')==c['table_sha256'],'changed table')
    for name,h in s['generator_files']:
        require(digest(ROOT/name)==h,'changed certificate data '+name)
    require(digest(ROOT/'replay/compiled/libmubk.dylib')==c['kernel_sha256'],'changed recorded binary')
    rows=ROOT/'replay'/f'{family}.jsonl'
    require(rows.is_file(),'missing detailed receipts')
    for line in rows.open():
        require(json.loads(line)['replay_context']==c,'mixed replay contexts')
    print(f'{family}: recorded identity and complete-summary fields PASS')
print('Integrity only: this does not authenticate status bits or replace clean mathematical replay.')
