"""Clean full replay, then repair unresolved parents by depth-three subdivision.

Usage: python3 chk/full_replay.py NEW_OUTPUT_DIRECTORY [workers]
No supplied replay receipt is imported. Exit zero requires both complete covers.
"""
import concurrent.futures
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent

def run(cmd, log):
    with log.open('w') as f:
        return subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT).returncode

def main():
    out = Path(sys.argv[1]).resolve()
    out.mkdir(exist_ok=False)
    workers = int(sys.argv[2]) if len(sys.argv)>2 else 4
    if not 1 <= workers <= 16: raise ValueError('workers must be 1..16')
    for family, gens in [('fourier',['results/gen_fourier.jsonl']),
                         ('szollosi',['results/gen_vm.jsonl','results/gen_vm_tail.jsonl','results/gen_mac.jsonl'])]:
        receipt=out/(family+'.jsonl')
        cmd=[sys.executable,'chk/run_check.py',family,'gen/table17.npy',str(receipt),str(workers),*gens]
        code=run(cmd+['--clean'],out/(family+'-base.log'))
        if code not in (0,2): raise RuntimeError(f'{family} base failed; see log')
        summary=json.loads(receipt.with_name(family+'_summary.json').read_text())
        if summary['cover_complete']:
            print(f'{family}: clean full replay PASS',flush=True);continue
        failed=[]
        for line in receipt.open():
            for box in json.loads(line)['boxes']:
                if box['ok'] is not True: failed.append([box['lev'],box['idx']])
        if not failed: raise RuntimeError('incomplete cover without repairable parents')
        print(f'{family}: {len(failed)} parents require subdivision',flush=True)
        jobs=[]
        for i,box in enumerate(failed):
            boxes=out/f'{family}-split-{i}.json'; boxes.write_text(json.dumps([box]))
            result=out/f'{family}-split-{i}.jsonl'
            scmd=[sys.executable,'chk/verify_split.py',family,'gen/table17.npy',str(result),*gens,'--boxes',str(boxes),'--depth','3']
            jobs.append((scmd,out/f'{family}-split-{i}.log',result))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures=[pool.submit(run,cmd,log) for cmd,log,result in jobs]
            for future in futures:
                if future.result()!=0: raise RuntimeError('subdivision subprocess failed')
        with receipt.open('ab') as f:
            for cmd,log,result in jobs: f.write(result.read_bytes())
        final_cmd=[sys.executable,'chk/run_check.py',family,'gen/table17.npy',str(receipt),str(workers),*gens]
        if run(final_cmd,out/(family+'-final.log')):
            raise RuntimeError(f'{family}: full cover not verified')
        print(f'{family}: full replay plus bound subdivisions PASS',flush=True)

if __name__=='__main__':main()
