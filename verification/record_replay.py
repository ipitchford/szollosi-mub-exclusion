"""Import a completed local replay and describe it, without authenticating it.

Usage: python3 verification/record_replay.py /path/to/completed-new-replay
Run only after the clean replay driver exits successfully. This is packaging,
not a replacement for that execution or a proof of honest receipt contents.
"""
import collections
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
source = Path(sys.argv[1]).resolve()
target = ROOT / 'replay'
if target.exists():
    raise RuntimeError('Refusing to replace an existing replay directory')
summaries = {}
for family in ('fourier', 'szollosi'):
    s = json.loads((source / f'{family}_summary.json').read_text())
    if not (s['cover_complete'] and s['verified_fail'] == 0
            and s['not_yet_checked'] == 0 and s['leaves_without_record'] == 0
            and s['verified_ok'] == s['boxes_needed']):
        raise RuntimeError(f'Incomplete {family} replay')
    summaries[family] = s
binary = ROOT / 'chk/libmubk.dylib'
digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
if any(s['replay_context']['kernel_sha256'] != digest(binary) for s in summaries.values()):
    raise RuntimeError('Current binary is not the recorded replay binary')
shutil.copytree(source, target)
(target / 'compiled').mkdir()
shutil.copy2(binary, target / 'compiled/libmubk.dylib')
subprocess.run([sys.executable, str(ROOT / 'verification/check_receipt_integrity.py')], check=True)
levels = collections.Counter()
for p in sorted(target.glob('szollosi-split-*.jsonl')):
    for line in p.read_text().splitlines():
        for b in json.loads(line)['boxes']:
            if b['ok'] is not True:
                raise RuntimeError('Failed subdivision receipt')
            levels[b['split_detail']['split_level']] += 1
s = summaries['szollosi']
c = s['replay_context']
now = datetime.datetime.now(datetime.timezone.utc).isoformat()
lines = [
    '# Local clean replay receipt', '',
    f'Recorded: {now}. Producer-side execution on one local machine.', '',
    'The revised clean driver was run into a new output directory. No supplied',
    'historical success record was imported. Required boundary subdivisions were',
    'computed afresh and then combined with this same run for final cover checking.', '',
    '```sh', 'sh -e chk/build.sh',
    'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 chk/full_replay.py NEW_OUTPUT_DIRECTORY 4',
    '```', '', '| Family | Required boxes | Verified | Cover |', '|---|---:|---:|---|'
]
for family, item in summaries.items():
    lines.append(f"| {family} | {item['boxes_needed']} | {item['verified_ok']} | complete |")
lines += ['', f'Szöllősi boundary parents replayed by subdivision: {sum(levels.values())}.',
          'Successful subdivision levels: ' + ', '.join(f'{n} parents at level {level}' for level, n in sorted(levels.items())) + '.',
          'Both final summaries have zero failed, unchecked or missing required leaves.', '',
          '## Bound evidence', '',
          '- `replay/fourier_summary.json` and `replay/szollosi_summary.json`: full source, table, certificate and environment identity.',
          '- `replay/*.jsonl`: per-root/per-box input-bound results; historical failed parent attempts remain visible before successful subdivisions.',
          '- `replay/*-base.log`, `replay/*-split-*.log`, and the final Szöllősi log: observed computation records.',
          '- `replay/compiled/libmubk.dylib`: exact macOS binary used for this run, not a portable executable.', '',
          f"Kernel SHA-256: `{c['kernel_sha256']}`.",
          f"Table SHA-256: `{c['table_sha256']}`.",
          f"Runtime: {c['python']}; NumPy {c['numpy']}; mpmath {c['mpmath']}; {c['platform']}.", '',
          '## Assurance ceiling', '',
          'This is local replay, not external independent reproduction, human peer',
          'review or end-to-end formal verification. Written ownership, symmetry and',
          'encoding arguments, the interval kernel, compiler and hardware remain',
          'separate trusted layers. Hashes bind bytes; they do not authenticate an',
          'arbitrary fabricated pass record. `check_receipt_integrity.py` checks',
          'recorded identity and summary completeness only; to reproduce the actual',
          'mathematical computation run the clean driver above.', '',
          'The manuscript distinguishes the new Szöllősi-family candidate from the',
          'Fourier theorem already proved in 2009. The unrestricted dimension-six',
          'three-basis conclusion remains conditional on MMVW Conjecture 2.', '']
(ROOT / 'REPLAY_RECEIPT.md').write_text('\n'.join(lines))
print('Imported completed replay and generated REPLAY_RECEIPT.md; no external independence asserted.')
