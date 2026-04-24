#!/bin/bash

logK_list2=(-0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6)
a_list2=(0.1 0.1 0.2 0.3 0.4 0.6 1.0 1.0 1.0 1.0 1.0 1.0 1.0)

######################################
### System A gamma/(M*omega_c)=32 ###
######################################
mkdir _sys_A_gam_32
cd _sys_A_gam_32

### Adiabatic BO rates ###
mkdir adiabatic
cd adiabatic
sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=A --meth=adi --fix=s --gam=32.0 --logK=-7.00 #g" ../../submit_template.slm
sbatch ../../submit_template.slm
for logK in "${logK_list2[@]}"; do
  sed -i "s#python.*#python ../../1_kcrpmd_tst.py --sys=A --meth=adi --fix=s --gam=32.0 --logK=$logK #g" ../../submit_template.slm
  sbatch ../../submit_template.slm
done
cd ../
### Adiabatic BO rates ###

cd ../
######################################
### System A gamma/(M*omega_c)=32 ###
######################################
