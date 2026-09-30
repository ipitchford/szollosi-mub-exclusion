#!/bin/bash
# wait for the generator (driver PID 73039 and its workers), then replay the remaining certificates of this machine
# and the tail certificates of machine A (6 processes); resumes from chk_szollosi_mac.jsonl
while kill -0 73039 2>/dev/null; do sleep 20; done
sleep 10
cd /Users/admin/mub6_szollosi/v2/chkmac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
/opt/homebrew/bin/python3.14 -u chk/run_check.py szollosi gen/table17.npy /Users/admin/mub6_szollosi/v2/chk_results/chk_szollosi_mac.jsonl 6 /Users/admin/mub6_szollosi/v2/results/gen_mac.jsonl /Users/admin/mub6_szollosi/v2/results/gen_vm_tail_copy.jsonl > /Users/admin/mub6_szollosi/v2/logs/chk_szollosi_mac_2.log 2>&1
