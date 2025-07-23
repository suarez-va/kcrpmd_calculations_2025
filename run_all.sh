#!/bin/bash
K0_list=(0.0000095 0.0001710 0.0005320 0.0009500 0.0016911 0.0030021 0.0094909)

#python kcrpmd_thermalization.py --sys 1 --method 1 --fix s --K0 0.00285 --leps 0.015 --hw 0

# System A BO rate at K=0
sed -i "s/python.*/python kcrpmd_thermalization.py --sys 1 --method 1 --fix s --K0 1e-10 --leps 0.015 --hw 0 /g" submit_template.slm
sbatch submit_template.slm 

# System A BO rates
for K0 in "${K0_list[@]}"; do
  sed -i "s/python.*/python kcrpmd_thermalization.py --sys 1 --method 1 --fix s --K0 $K0 --leps 0.015 --hw 0 /g" submit_template.slm
  sbatch submit_template.slm 
done

# System A KC-RPMD rates
for method in 2 3 ; do
  for fix in s y ; do
    for K0 in "${K0_list[@]}"; do
      sed -i "s/python.*/python kcrpmd_thermalization.py --sys 1 --method $method --fix $fix --K0 $K0 --leps 0.015 --hw 0 /g" submit_template.slm
      sbatch submit_template.slm 
    done
  done
done

