#!/bin/bash

ntraj=1000
check_every=33
submitted_since_check=0
njobs=0

for dir in _sys_3*/; do
  cd "$dir/libra_data" || continue
  pwd

  for ((itraj=1; itraj<=ntraj; itraj++)); do
    while true; do
      [[ -d _itraj_$itraj ]] && break

      # Only check the queue every N submissions
      if (( submitted_since_check % check_every == 0 )); then
        njobs=$(squeue -u "$USER" | tail -n +2 | wc -l)
        submitted_since_check=0
      fi

      if (( njobs < 399 )); then
        sed -i "s|python.*|python ../../3_kcrpmd_dynamics.py --itraj=$itraj|g" ../../submit_template_sys3.slm
        echo "$itraj"
        sbatch ../../submit_template_sys3.slm

        ((submitted_since_check++))
        break
      else
        echo "hold on one sec..."
        sleep 10
      fi
    done
  done

  cd ../../
done

