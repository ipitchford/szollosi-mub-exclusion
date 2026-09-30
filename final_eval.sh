#!/bin/bash
# final_eval.sh -- merge the checker outputs, run the final cover evaluation, and compute every statistic of the paper.
set -e
echo 'Historical v2 driver disabled: use chk/full_replay.py with a new output directory for revision 3.' >&2
exit 2
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# 1. merged Szollosi checker output: machine B's replays first, then machine A's parameter-split verifications
cat chk_results/chk_szollosi_mac.jsonl chk_results/chk_szollosi_vmrec_mac.jsonl chk_results/chk_szollosi_split.jsonl chk_results/chk_szollosi_split_b.jsonl chk_results/chk_szollosi_split_mac2_?.jsonl > chk_results/chk_szollosi_all.jsonl
# 2. final cover evaluation (every box is already verified; the run only re-verifies the table and evaluates the cover tree)
python3 -u chk/run_check.py szollosi gen/table17.npy chk_results/chk_szollosi_all.jsonl 2 results/gen_vm.jsonl results/gen_vm_tail.jsonl results/gen_mac.jsonl > logs/chk_szollosi_final.log 2>&1
cp chk_results/chk_szollosi_all_summary.json chk_results/chk_szollosi_summary.json
# 3. statistics
python3 analysis/stats.py szollosi analysis/stats_szollosi.json --gen results/gen_vm.jsonl results/gen_vm_tail.jsonl results/gen_mac.jsonl --chk chk_results/chk_szollosi_all.jsonl > /dev/null
python3 analysis/stats.py fourier analysis/stats_fourier.json --gen results/gen_fourier.jsonl --chk chk_results/chk_fourier.jsonl > /dev/null
python3 analysis/special_boxes.py analysis/special_szollosi.json chk_results/chk_szollosi_all.jsonl --gen results/gen_vm.jsonl results/gen_vm_tail.jsonl results/gen_mac.jsonl
python3 analysis/empty_boxes.py szollosi analysis/empty_szollosi.json chk_results/chk_szollosi_all.jsonl
python3 paper/make_values.py > /dev/null
python3 paper/make_values2.py > /dev/null
echo "final evaluation done"
