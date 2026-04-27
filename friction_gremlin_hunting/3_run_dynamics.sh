#!/bin/bash

ntraj=9999
check_every=20
submitted_since_check=0
njobs=0

for dir1 in _sys_A_gam_32/; do
  cd $dir1

  if [ -d "adiabatic" ]; then
    cd adiabatic
    for dir2 in _fix_*/; do
      cd "$dir2"
      mkdir libra_data
      cd libra_data
      #cd "$dir2/libra_data" || continue
      pwd
      for ((itraj=1; itraj<=ntraj; itraj++)); do
        while true; do
          [ -f "_itraj_$itraj/mem_data.hdf" ] && break 
          # Only check the queue every N submissions
          if (( submitted_since_check % check_every == 0 )); then
            njobs=$(squeue -u "$USER" | tail -n +2 | wc -l)
            submitted_since_check=0
          fi
          if (( njobs < 480 )); then
            sed -i "s|python.*|python ../../../../3_kcrpmd_dynamics.py --itraj=$itraj|g" ../../../../submit_template.slm
            sbatch ../../../../submit_template.slm
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
done

