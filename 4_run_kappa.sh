#!/bin/bash

for dir in _sys_1_gam_32*/; do
  cd $dir
  pwd
  sed -i "s|python.*|python ../4_kcrpmd_kappa.py|g" ../submit_template.slm
  sbatch ../submit_template.slm
  cd ../
done

