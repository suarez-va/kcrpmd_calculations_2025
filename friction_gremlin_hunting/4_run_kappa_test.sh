#!/bin/bash

for dir1 in _sys_*/; do
  cd $dir1
  pwd

  if [ -d "adiabatic" ]; then
    cd adiabatic
    cd "_fix_s_logK_-0.10"
    pwd
    sed -i "s|python.*|python ../../../4_kcrpmd_kappa.py |g" ../../../submit_template_test.slm
    sbatch ../../../submit_template_test.slm
    cd ../
    cd ../
  fi



  cd ../
done

