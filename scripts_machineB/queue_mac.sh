#!/bin/bash
# wait for the first checker run, then replay the VM's certificates on this machine (2 processes)
while kill -0 48397 2>/dev/null; do sleep 30; done
cd /Users/admin/mub6_szollosi/v2/chkmac
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
/opt/homebrew/bin/python3.14 -u chk/run_check.py szollosi gen/table17.npy /Users/admin/mub6_szollosi/v2/chk_results/chk_szollosi_vmrec_mac.jsonl 2 /Users/admin/mub6_szollosi/v2/results/gen_vm_copy.jsonl > /Users/admin/mub6_szollosi/v2/logs/chk_szollosi_vmrec_mac.log 2>&1
