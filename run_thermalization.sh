#!/bin/bash

#python kcrpmd_thermalization.py --sys=1 --method=1 --fix=s --K0=0.00285 --leps=0.015 --hw=0
#K0_list=(9.50e-06 9.50e-05 3.00e-04 5.34e-04 9.50e-04 1.69e-03 3.00e-03 5.34e-03)
K0_list=(9.50e-06 9.50e-05 3.00e-04 5.34e-04 9.50e-04 1.69e-03 3.00e-03 4.50e-03 6.00e-03)
leps_list=(-1.24e-02 -1.33e-02 -1.43e-02 -1.52e-02 -1.62e-02)

# System A BO rate at K=0
sed -i "s/python.*/python kcrpmd_thermalization.py --sys=1 --method=1 --fix=s --K0=1e-10 --leps=-1.43e-02 --hw=0 /g" submit_template.slm
sbatch submit_template.slm 

# System A BO rates
for K0 in "${K0_list[@]}"; do
  sed -i "s/python.*/python kcrpmd_thermalization.py --sys=1 --method=1 --fix=s --K0=$K0 --leps=-1.43e-02 --hw=0 /g" submit_template.slm
  sbatch submit_template.slm 
done

# System A KC-RPMD rates
for method in 2 3 ; do
  for fix in s y ; do
    for K0 in "${K0_list[@]}"; do
      sed -i "s/python.*/python kcrpmd_thermalization.py --sys=1 --method=$method --fix=$fix --K0=$K0 --leps=-1.43e-02 --hw=0 /g" submit_template.slm
      sbatch submit_template.slm 
    done
  done
done

# System B BO rate at K=0
sed -i "s/python.*/python kcrpmd_thermalization.py --sys=2 --method=1 --fix=s --K0=1e-10 --leps=-1.43e-02 --hw=0 /g" submit_template.slm
sbatch submit_template.slm 

# System B BO rates
for K0 in "${K0_list[@]}"; do
  sed -i "s/python.*/python kcrpmd_thermalization.py --sys=2 --method=1 --fix=s --K0=$K0 --leps=-1.43e-02 --hw=0 /g" submit_template.slm
  sbatch submit_template.slm 
done

# System B KC-RPMD rates
for method in 2 3 ; do
  for fix in s y ; do
    for K0 in "${K0_list[@]}"; do
      sed -i "s/python.*/python kcrpmd_thermalization.py --sys=2 --method=$method --fix=$fix --K0=$K0 --leps=-1.43e-02 --hw=0 /g" submit_template.slm
      sbatch submit_template.slm 
    done
  done
done

# System C BO rate at K=0 w/ adiabatic and non-adiabatic hard wall rates
for leps in "${leps_list[@]}"; do
  for hw in 1 -1 ; do
    sed -i "s/python.*/python kcrpmd_thermalization.py --sys=3 --method=1 --fix=s --K0=1e-10 --leps=$leps --hw=$hw /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System C BO adiabatic and non-adiabatic hard wall rates
for leps in "${leps_list[@]}"; do
  for hw in 1 -1 ; do
    sed -i "s/python.*/python kcrpmd_thermalization.py --sys=3 --method=1 --fix=s --K0=2.85e-03 --leps=$leps --hw=$hw /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System C KC-RPMD adiabatic hard wall rates
for method in 2 3 ; do
  for leps in "${leps_list[@]}"; do
    sed -i "s/python.*/python kcrpmd_thermalization.py --sys=3 --method=$method --fix=s --K0=2.85e-03 --leps=$leps --hw=1 /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System C KC-RPMD nonadiabatic hard wall rates
for method in 2 3 ; do
  for leps in "${leps_list[@]}"; do
    sed -i "s/python.*/python kcrpmd_thermalization.py --sys=3 --method=$method --fix=y --K0=2.85e-03 --leps=$leps --hw=-1 /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

