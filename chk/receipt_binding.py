"""Content binding for local replay receipts; hashes are not authentication."""
import hashlib
import json
import os
import platform
import sys
from pathlib import Path
import numpy as np
import mpmath

HERE = Path(__file__).resolve().parent

def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()

def context(fam, table):
    if os.environ.get('MUB6_KERNEL'):
        raise RuntimeError('MUB6_KERNEL override forbidden in proof replay')
    kernel = HERE / ('libmubk.dylib' if sys.platform == 'darwin' else 'libmubk.so')
    return {
        'schema': 'mub6-replay/v3', 'family': fam.name,
        'family_data': {'E': fam.E.tolist(), 'coefficients': fam.coef_rect.tolist(),
                        'CE': fam.CE.tolist(), 'constraint_coefficients': fam.ct_rect.tolist()},
        'table_sha256': sha256(table),
        'code_sha256': {n: sha256(HERE / n) for n in
                       ('ck.py', 'kernel.c', 'run_check.py', 'verify_split.py', 'receipt_binding.py', 'full_replay.py', 'build.sh')},
        'kernel_sha256': sha256(kernel), 'python': sys.version,
        'numpy': np.__version__, 'mpmath': mpmath.__version__, 'platform': platform.platform(),
    }

def binding(ctx, key, hulls):
    lev, idx = key
    return digest({'context': ctx, 'parameter_box': [lev, list(idx)], 'hulls': hulls})

def read_cache(path, ctx, recs):
    done = {}
    if not Path(path).exists():
        return done
    with open(path) as f:
        for line in f:
            row = json.loads(line)
            if row.get('replay_context') != ctx:
                raise RuntimeError('unbound or mismatched replay context; use a new output path')
            for box in row.get('boxes', []):
                key = (box['lev'], tuple(box['idx']))
                if key not in recs or box.get('input_sha256') != binding(ctx, key, recs[key]):
                    raise RuntimeError('cached certificate inputs changed; use a new output path')
                if type(box.get('ok')) is not bool:
                    raise RuntimeError('invalid replay status')
                done[key] = box['ok']
    return done
