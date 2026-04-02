#!/bin/bash

logK_list1=(-2.00 -1.00 -0.50 -0.25 0.00 0.25 0.50 0.66 0.82 0.98)
a_list1=(0.1 0.1 0.1 0.3 1.0 1.0 1.0 1.0 1.0 1.0)

logK_list2=(-0.90 -0.75 -0.60 -0.45 -0.30 -0.15 0.00 0.15 0.30 0.45 0.60)
a_list2=(0.1 0.1 0.1 0.2 0.3 0.5 1.0 1.0 1.0 1.0 1.0)

leps_list=(0.00 4.00 8.00 12.00 16.00)

####################################
### System A gamma/(M*omega_c)=1 ###
####################################
mkdir _sys_A_gam_1
cd _sys_A_gam_1

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=1.0 --logK=-7.00 |g" ../../submit_template.slm
sbatch ../../submit_template.slm
for logK in "${logK_list2[@]}"; do
  sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=1.0 --logK=$logK |g" ../../submit_template.slm
  sbatch ../../submit_template.slm
done
cd ../
### Adiabatic BO rates ###

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list2[@]}"; do
    sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=1 --meth=3 --fix=$fix --a=0.1 --gam=1.0 --logK=$logK |g" ../../submit_template.slm
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
### System A gamma/(M*omega_c)=32 ###
######################################
mkdir _sys_A_gam_32
cd _sys_A_gam_32

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=32.0 --logK=-7.00 |g" ../../submit_template.slm
sbatch ../../submit_template.slm
for logK in "${logK_list2[@]}"; do
  sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=1 --meth=1 --fix=s --gam=32.0 --logK=$logK |g" ../../submit_template.slm
  sbatch ../../submit_template.slm
done
cd ../
### Adiabatic BO rates ###

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list2[@]}"; do
    sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=1 --meth=3 --fix=$fix --a=0.1 --gam=32.0 --logK=$logK |g" ../../submit_template.slm
    sbatch ../../submit_template.slm
  done
done
cd ../
### New KC-RPMD rates ###

cd ../
######################################
### System A gamma/(M*omega_c)=32 ###
######################################

