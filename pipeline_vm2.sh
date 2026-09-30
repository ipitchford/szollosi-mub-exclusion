#!/bin/bash
# VM pipeline, resumed after a restart of the cloud container (the generator resumes from its results file):
#  1. Szollosi generator on the last 250 roots of machine B's list (tail-first), 2. Fourier generator, 3. Fourier checker
cd /home/claude/mub6/v2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "$(date) container restarted; pipeline resumed" >> logs/pipeline_vm.log
cd gen && MUB6_BASE=/home/claude/mub6/v2 MUB6_RESULTS=gen_vm_tail.jsonl MUB6_REVERSE=1 python3 -u driver.py szollosi 2 0,1,2,3,4,5/7 250 > ../logs/gen_vm_tail2.log 2>&1; cd ..
echo "$(date) szollosi generator (tail of the other list) finished" >> logs/pipeline_vm.log
cd gen && MUB6_BASE=/home/claude/mub6/v2 MUB6_RESULTS=gen_fourier.jsonl python3 -u driver.py fourier 2 > ../logs/gen_fourier.log 2>&1; cd ..
echo "$(date) fourier generator finished" >> logs/pipeline_vm.log
python3 -u chk/run_check.py fourier gen/table17.npy chk_results/chk_fourier.jsonl 2 results/gen_fourier.jsonl > logs/chk_fourier.log 2>&1
echo "$(date) fourier checker finished" >> logs/pipeline_vm.log
