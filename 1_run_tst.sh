#!/bin/bash

logK_list1=(-2.00 -1.00 -0.50 -0.25 0.00 0.25 0.50 0.66 0.82 0.98)
a_list1=(0.1 0.1 0.1 0.3 1.0 1.0 1.0 1.0 1.0 1.0)

logK_list2=(-0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6)
a_list2=(0.1 0.1 0.2 0.3 0.4 0.6 1.0 1.0 1.0 1.0 1.0 1.0 1.0)

leps_list=(0.00 4.00 8.00 12.00 16.00)

####################################
### System A gamma/(M*omega_c)=1 ###
####################################
mkdir _sys_A_gam_1
cd _sys_A_gam_1

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=1.0 --logK=-7.00 #g" ../../submit_template.slm
#sbatch ../../submit_template.slm
for logK in "${logK_list1[@]}"; do
  sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=1.0 --logK=$logK #g" ../../submit_template.slm
  sbatch ../../submit_template.slm
done
cd ../
### Adiabatic BO rates ###

### Original KC-RPMD rates ###
mkdir kcrpmd_ori
cd kcrpmd_ori
for fix in s y; do
  for i in "${!logK_list1[@]}"; do
    logK=${logK_list1[i]}
    a=${a_list1[i]}
    sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=2 --fix=$fix --a=$a --gam=1.0 --logK=$logK #g" ../../submit_template.slm
    sbatch ../../submit_template.slm
    if [[ "$a" != "0.1" ]]; then
      sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=2 --fix=$fix --a=0.1 --gam=1.0 --logK=$logK #g" ../../submit_template.slm
      sbatch ../../submit_template.slm
    fi
  done
done
cd ../
### Original KC-RPMD rates ###

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list1[@]}"; do
    sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=3 --fix=$fix --a=0.1 --gam=1.0 --logK=$logK #g" ../../submit_template.slm
    sbatch ../../submit_template.slm
  done
done
cd ../
### New KC-RPMD rates ###

cd ../
####################################
### System A gamma/(M*omega_c)=1 ###
####################################

exit 0

mkdir _sys_A_gam_250
mkdir _sys_B
mkdir _sys_C

# System A high friction BO rate at K=0
sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --gam=32 --method=1 --fix=s --K0=1e-10 /g" submit_template.slm
sbatch submit_template.slm

# System A high friction BO rates
for K0 in "${K0_list[@]}"; do
  sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --gam=32 --method=1 --fix=s --K0=$K0 /g" submit_template.slm
  sbatch submit_template.slm
done

# System A high friction original KC-RPMD rates
for fix in s y; do
  for i in "${!K0_list[@]}"; do
    K0=${K0_list[i]}
    a=${a_list[i]}
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --gam=32 --method=2 --fix=$fix --a=$a --K0=$K0 /g" submit_template.slm
    sbatch submit_template.slm
  done
done

# System A high friction original KC-RPMD rates set to a=0.1
for fix in s y; do
  for i in "${!K0_list[@]}"; do
    K0=${K0_list[i]}
    if [[ $i -lt 3 ]]; then
      echo $K0
      continue
    fi
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --gam=32 --method=2 --fix=$fix --K0=$K0 /g" submit_template.slm
    sbatch submit_template.slm
  done
done

# System A high friction new KC-RPMD rates
for fix in s y; do
  for K0 in "${K0_list[@]}"; do
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=1 --gam=32 --method=3 --fix=$fix --K0=$K0 /g" submit_template.slm
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
for fix in s y; do
  for K0 in "${K0_list[@]}"; do
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=2 --method=3 --fix=$fix --K0=$K0 /g" submit_template.slm
    sbatch submit_template.slm
  done
done

# System C BO rate at K=0 w/ adiabatic and non-adiabatic hard wall rates
for leps in "${leps_list[@]}"; do
  for hw in 1 -1; do
    sed -i "s/python.*/python 1_kcrpmd_tst.py --sys=3 --method=1 --fix=s --K0=1e-10 --leps=$leps --hw=$hw /g" submit_template.slm
    sbatch submit_template.slm
  done
done

# System C BO adiabatic and non-adiabatic hard wall rates
for leps in "${leps_list[@]}"; do
  for hw in 1 -1; do
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
