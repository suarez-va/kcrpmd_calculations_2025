#!/bin/bash

ntraj=1000

for dir in _sys_*/; do
  cd $dir
  pwd
  for ((itraj=1; itraj<=ntraj; itraj++)); do
    while true; do
      [[ -d _itraj_$itraj ]] && break
      njobs=$(squeue -u "$USER" | tail -n +2 | wc -l)
      if [ "$njobs" -lt 499 ]; then
        sed -i "s|python.*|python ../kcrpmd_transmission.py --itraj=$itraj|g" ../submit_template.slm
        echo "$itraj"
        sbatch ../submit_template.slm
        break
      else
        echo "hold on one sec..."
        sleep 60
      fi
    done
  done
  cd ../
done

