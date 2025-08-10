#!/bin/bash

for dir in _sys_2_*/; do
  cd $dir
  pwd
  sed -i "s|python.*|python ../kcrpmd_kappacalc.py|g" ../submit_template.slm
  sbatch ../submit_template.slm
  cd ../
done

