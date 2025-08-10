#!/bin/bash

#python 1_kcrpmd_tst.py --sys=1 --method=1 --fix=s --K0=0.00285 --leps=0.015 --hw=0
a_list=(0.1 0.1 0.1 0.3 1.0 1.0 1.0 1.0 1.0 1.0)
K0_list=(9.55e-06 9.55e-05 3.02e-04 5.37e-04 9.55e-04 1.70e-03 3.02e-03 4.47e-03 6.46e-03 9.33e-03)
leps_list=(-1.24e-02 -1.33e-02 -1.43e-02 -1.52e-02 -1.62e-02)

# System A BO rate at K=0
sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --method=1 --fix=s --K0=1e-10 /g" submit_template.slm
sbatch submit_template.slm 

# System A BO rates
for K0 in "${K0_list[@]}"; do
  sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --method=1 --fix=s --K0=$K0 /g" submit_template.slm
  sbatch submit_template.slm 
done

# System A original KC-RPMD rates
for fix in s y ; do
  for i in "${!K0_list[@]}"; do
    K0=${K0_list[i]}
    a=${a_list[i]}
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --method=2 --fix=$fix --a=$a --K0=$K0 /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System A original KC-RPMD rates set to a=0.1
for fix in s y ; do
  for i in "${!K0_list[@]}"; do
    K0=${K0_list[i]}
    if [[ $i -lt 3 ]]; then
      echo $K0
      continue
    fi
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --method=2 --fix=$fix --K0=$K0 /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System A new KC-RPMD rates
for fix in s y ; do
  for K0 in "${K0_list[@]}"; do
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --method=3 --fix=$fix --K0=$K0 /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System B BO rate at K=0
sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=2 --method=1 --fix=s --K0=1e-10 /g" submit_template.slm
sbatch submit_template.slm 

# System B BO rates
for K0 in "${K0_list[@]}"; do
  sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=2 --method=1 --fix=s --K0=$K0 /g" submit_template.slm
  sbatch submit_template.slm 
done

# System B new KC-RPMD rates
for fix in s y ; do
  for K0 in "${K0_list[@]}"; do
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=2 --method=3 --fix=$fix --K0=$K0 /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System C BO rate at K=0 w/ adiabatic and non-adiabatic hard wall rates
for leps in "${leps_list[@]}"; do
  for hw in 1 -1 ; do
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=3 --method=1 --fix=s --K0=1e-10 --leps=$leps --hw=$hw /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System C BO adiabatic and non-adiabatic hard wall rates
for leps in "${leps_list[@]}"; do
  for hw in 1 -1 ; do
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=3 --method=1 --fix=s --leps=$leps --hw=$hw /g" submit_template.slm
    sbatch submit_template.slm 
  done
done

# System C KC-RPMD adiabatic hard wall rates
for leps in "${leps_list[@]}"; do
  sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=3 --method=3 --fix=s --leps=$leps --hw=1 /g" submit_template.slm
  sbatch submit_template.slm 
done

# System C KC-RPMD nonadiabatic hard wall rates
for leps in "${leps_list[@]}"; do
  sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=3 --method=3 --fix=y --leps=$leps --hw=-1 /g" submit_template.slm
  sbatch submit_template.slm 
done

