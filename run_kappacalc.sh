#!/bin/bash

for dir in _sys_1_method_1*/; do
  cd $dir
  pwd
  sed -i "s|python.*|python ../kcrpmd_kappacalc.py|g" ../submit_template.slm
  sbatch ../submit_template.slm
  cd ../
done

