#!/bin/bash
# sequential VM pipeline (see logs/pipeline_vm.log):
#  1. wait for the VM share of the Szollosi generator (roots i = 6 mod 7)
#  2. Szollosi generator on the last 250 roots of the other machine's list (i != 6 mod 7), processed tail-first
#  3. Fourier generator, 4. Fourier checker, 5. checker replay of the VM's Szollosi certificates
cd /home/claude/mub6/v2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while kill -0 594 2>/dev/null; do sleep 30; done
echo "$(date) szollosi generator (VM share) finished" >> logs/pipeline_vm.log
cd gen && MUB6_BASE=/home/claude/mub6/v2 MUB6_RESULTS=gen_vm_tail.jsonl MUB6_REVERSE=1 python3 -u driver.py szollosi 2 0,1,2,3,4,5/7 250 > ../logs/gen_vm_tail.log 2>&1; cd ..
echo "$(date) szollosi generator (tail of the other list) finished" >> logs/pipeline_vm.log
cd gen && MUB6_BASE=/home/claude/mub6/v2 MUB6_RESULTS=gen_fourier.jsonl python3 -u driver.py fourier 2 > ../logs/gen_fourier.log 2>&1; cd ..
echo "$(date) fourier generator finished" >> logs/pipeline_vm.log
python3 -u chk/run_check.py fourier gen/table17.npy chk_results/chk_fourier.jsonl 2 results/gen_fourier.jsonl > logs/chk_fourier.log 2>&1
echo "$(date) fourier checker finished" >> logs/pipeline_vm.log
python3 -u chk/run_check.py szollosi gen/table17.npy chk_results/chk_szollosi_vm.jsonl 2 results/gen_vm.jsonl results/gen_vm_tail.jsonl > logs/chk_szollosi_vm.log 2>&1
echo "$(date) szollosi checker (VM records) finished" >> logs/pipeline_vm.log
