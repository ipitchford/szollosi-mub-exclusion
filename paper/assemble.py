"""Assemble the revision templates without importing historical replay results."""
import json
from pathlib import Path
import re

root = Path(__file__).resolve().parent
values = json.loads((root / 'values.json').read_text())
parts = [str(i) for i in range(1, 12)] + ['A', 'B', 'C', 'S']
text = '\n\n'.join((root / 'sections' / f'part{i}.md').read_text() for i in parts)
for key, value in values.items():
    text = text.replace('@@' + key + '@@', str(value))
missing = re.findall(r'@@[A-Z0-9_]+@@', text)
if missing:
    raise ValueError(f'Unfilled placeholders: {missing}')
(root / 'paper.md').write_text(text)
