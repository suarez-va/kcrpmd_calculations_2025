#!/bin/bash

ntraj=1000

for dir in _sys_*/; do
  cd $dir
  pwd
  for ((itraj=1; itraj<=ntraj; itraj++)); do
    while true; do
      [[ -d _itraj_$itraj ]] && break
      njobs=$(squeue -u "$USER" | tail -n +2 | wc -l)
      if [ "$njobs" -lt 400 ]; then
        sed -i "s|python.*|python ../kcrpmd_transmission.py --itraj=$itraj --nsteps=100000|g" ../submit_template.slm
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

#[[ $(basename "$PWD") == _sys_* ]] && echo "$PWD" || { echo "not a thermalization directory"; exit; }
#
#ntraj=1000
#for ((itraj=1; itraj<=ntraj; itraj++)); do
#  while true; do
#    [[ -d _itraj_$itraj ]] && break
#    njobs=$(squeue -u "$USER" | tail -n +2 | wc -l)
#    if [ "$njobs" -lt 200 ]; then
#      sed -i "s|python.*|python ../kcrpmd_transmission.py --itraj=$itraj --nsteps=100000|g" ../submit_template.slm
#      echo "$itraj"
#      sbatch ../submit_template.slm
#      break
#    else
#      echo "$njobs"
#      sleep 60
#    fi
#  done
#done
