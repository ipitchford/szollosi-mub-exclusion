"""fill.py -- substitute the @@...@@ placeholders of the paper and the response from the analysis outputs.
Usage: python3 fill.py VALUES.json IN.md OUT.md  (fails if a placeholder remains)"""
import sys, json, re
vals = json.load(open(sys.argv[1]))
s = open(sys.argv[2]).read()
for k, v in vals.items():
    s = s.replace('@@' + k + '@@', str(v))
left = sorted(set(re.findall(r'@@[A-Z0-9_]+@@', s)))
if left:
    print('UNFILLED:', left); sys.exit(1)
open(sys.argv[3], 'w').write(s)
print('filled', sys.argv[3])
