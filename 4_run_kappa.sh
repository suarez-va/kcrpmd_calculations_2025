#!/bin/bash

for dir in _sys_*/; do
  cd $dir
  pwd
  sed -i "s|python.*|python ../4_kcrpmd_kappa.py|g" ../submit_template.slm
  sbatch ../submit_template.slm
  cd ../
done

