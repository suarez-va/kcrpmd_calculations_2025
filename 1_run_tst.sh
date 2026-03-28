#!/bin/bash

logK_list1=(-2.00 -1.00 -0.50 -0.25 0.00 0.25 0.50 0.66 0.82 0.98)
a_list1=(0.1 0.1 0.1 0.3 1.0 1.0 1.0 1.0 1.0 1.0)

logK_list2=(-0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6)
a_list2=(0.1 0.1 0.2 0.3 0.4 0.6 1.0 1.0 1.0 1.0 1.0 1.0 1.0)

leps_list=(0.00 4.00 8.00 12.00 16.00)

################
### System A ###
################
mkdir _sys_A
cd _sys_A

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=1.0 --logK=-7.00 #g" ../../submit_template.slm
sbatch ../../submit_template.slm
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
################
### System A ###
################

####################################
### System A gamma/(M*omega_c)=1 ###
####################################
mkdir _sys_A_gam_1
cd _sys_A_gam_1

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=1.0 --logK=-7.00 #g" ../../submit_template.slm
sbatch ../../submit_template.slm
for logK in "${logK_list2[@]}"; do
  sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=1.0 --logK=$logK #g" ../../submit_template.slm
  sbatch ../../submit_template.slm
done
cd ../
### Adiabatic BO rates ###

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list2[@]}"; do
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

######################################
### System A gamma/(M*omega_c)=250 ###
######################################
mkdir _sys_A_gam_250
cd _sys_A_gam_250

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=250.0 --logK=-7.00 #g" ../../submit_template.slm
sbatch ../../submit_template.slm
for logK in "${logK_list2[@]}"; do
  sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=250.0 --logK=$logK #g" ../../submit_template.slm
  sbatch ../../submit_template.slm
done
cd ../
### Adiabatic BO rates ###

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list2[@]}"; do
    sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=1 --meth=3 --fix=$fix --a=0.1 --gam=250.0 --logK=$logK #g" ../../submit_template.slm
    sbatch ../../submit_template.slm
  done
done
cd ../
### New KC-RPMD rates ###

cd ../
######################################
### System A gamma/(M*omega_c)=250 ###
######################################

################
### System B ###
################
mkdir _sys_B
cd _sys_B

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=2 --meth=1 --fix=s --gam=1.0 --logK=-7.00 #g" ../../submit_template.slm
sbatch ../../submit_template.slm
for logK in "${logK_list1[@]}"; do
  sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=2 --meth=1 --fix=s --gam=1.0 --logK=$logK #g" ../../submit_template.slm
  sbatch ../../submit_template.slm
done
cd ../
### Adiabatic BO rates ###

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list1[@]}"; do
    sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=2 --meth=3 --fix=$fix --a=0.1 --gam=1.0 --logK=$logK #g" ../../submit_template.slm
    sbatch ../../submit_template.slm
  done
done
cd ../
### New KC-RPMD rates ###

cd ../
################
### System B ###
################

################
### System C ###
################
mkdir _sys_C
cd _sys_C

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
for leps in "${leps_list[@]}"; do
  for hw in 1 -1; do
    sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=3 --meth=1 --fix=s --gam=1.0 --logK=-7.00 --leps=$leps --hw=$hw #g" ../../submit_template.slm
    sbatch ../../submit_template.slm
    sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=3 --meth=1 --fix=s --gam=1.0 --logK=0.50 --leps=$leps --hw=$hw #g" ../../submit_template.slm
    sbatch ../../submit_template.slm
  done
done
cd ../
### Adiabatic BO rates ###

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for leps in "${leps_list[@]}"; do
  for hw in 1 -1; do
    if [[ "$hw" -= "1" ]]; then
      sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=3 --meth=3 --fix=s --gam=1.0 --logK=0.50 --leps=$leps --hw=$hw #g" ../../submit_template.slm
      sbatch ../../submit_template.slm
    fi
    if [[ "$hw" -= "-1" ]]; then
      sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=3 --meth=3 --fix=y --gam=1.0 --logK=0.50 --leps=$leps --hw=$hw #g" ../../submit_template.slm
      sbatch ../../submit_template.slm
    fi
  done
done
cd ../
### New KC-RPMD rates ###

cd ../
################
### System C ###
################

