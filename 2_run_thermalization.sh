#!/bin/bash

for dir in _sys_*/; do
  cd $dir
  pwd
  sed -i "s|python.*|python ../2_kcrpmd_thermalization.py|g" ../submit_template.slm
  sbatch ../submit_template.slm
  cd ../
done

