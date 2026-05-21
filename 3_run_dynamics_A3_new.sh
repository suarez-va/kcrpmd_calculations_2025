#!/bin/bash

ntraj=12250
check_every=5
submitted_since_check=0
njobs=0

cd _sys_A3

if [ -d "kcrpmd_new" ]; then
  cd kcrpmd_new
  for dir2 in _fix_*/; do
    cd "$dir2/libra_data" || continue
    pwd
    for ((itraj=1; itraj<=ntraj; itraj++)); do
      while true; do
        [ -f "_itraj_$itraj/mem_data.hdf" ] && break 
        echo $itraj
        # Only check the queue every N submissions
        if (( submitted_since_check % check_every == 0 )); then
          njobs=$(squeue -u "$USER" | tail -n +2 | wc -l)
          submitted_since_check=0
        fi
        if (( njobs < 489 )); then
          sed -i "s|python.*|python ../../../../3_kcrpmd_dynamics.py --itraj=$itraj --iskip=4|g" ../../../../submit_template_A3_new.slm
          sbatch ../../../../submit_template_A3_new.slm
          ((submitted_since_check++))
          break
        else
          echo "hold on one sec..."
          sleep 5
        fi
      done
    done
    cd ../../
  done
  cd ../
fi

cd ../

