"""Validate exact file coverage, not just listed hashes."""
import hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def excluded(p):
    return (any(x in {'.git','__pycache__'} for x in p.parts)
            or str(p) in {'chk/libmubk.dylib','chk/libmubk.so'}
            or p.suffix in {'.pyc','.aux','.log','.out'} and p.parts[0]=='paper')
expected={}
for line in (ROOT/'MANIFEST.sha256').read_text().splitlines():
    h,name=line.split('  ',1)
    if name in expected: raise RuntimeError('duplicate manifest entry')
    expected[name]=h
actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and not excluded(p.relative_to(ROOT)) and p.name not in {'MANIFEST.sha256','SHA256SUMS'}}
if set(expected)!=actual: raise RuntimeError(f'coverage mismatch: missing={actual-set(expected)}, extra={set(expected)-actual}')
for name,h in expected.items():
    if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=h: raise RuntimeError('hash mismatch: '+name)
if (ROOT/'SHA256SUMS').read_bytes()!=(ROOT/'MANIFEST.sha256').read_bytes(): raise RuntimeError('manifest aliases differ')
print(f'PASS: exact coverage and SHA-256 for {len(expected)} files')
