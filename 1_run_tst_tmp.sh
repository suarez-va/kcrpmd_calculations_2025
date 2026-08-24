#!/bin/bash

logK_list1=(-2.00 -1.00 -0.50 -0.25 0.00 0.25 0.50 0.66 0.82 0.98)
a_list1=(0.1 0.1 0.1 0.3 1.0 1.0 1.0 1.0 1.0 1.0)

############################
#### System A1 old gamma ###
############################
#mkdir _sys_A1_oldgamma
#cd _sys_A1_oldgamma
#
#### Original KC-RPMD rates ###
#mkdir kcrpmd_ori
#cd kcrpmd_ori
#for fix in s y; do
#  for i in "${!logK_list1[@]}"; do
#    logK=${logK_list1[i]}
#    a=${a_list1[i]}
#    sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=A1 --meth=ori --fix=$fix --a=$a --logK=$logK |g" ../../submit_template.slm
#    sbatch ../../submit_template.slm
#    if [[ "$a" != "0.1" ]]; then
#      sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=A1 --meth=ori --fix=$fix --a=0.1 --logK=$logK |g" ../../submit_template.slm
#      sbatch ../../submit_template.slm
#    fi
#  done
#done
#cd ../
#### Original KC-RPMD rates ###
#
#### New KC-RPMD rates ###
#mkdir kcrpmd_new
#cd kcrpmd_new
#for fix in s y; do
#  for logK in "${logK_list1[@]}"; do
#    sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=A1 --meth=new --fix=$fix --a=0.1 --logK=$logK |g" ../../submit_template.slm
#    sbatch ../../submit_template.slm
#  done
#done
#cd ../
#### New KC-RPMD rates ###
#
#cd ../
############################
#### System A1 old gamma ###
############################

##############################
### System B a=1.0, c=1.0 ###
##############################
mkdir _sys_B_a1.0_c1.0
cd _sys_B_a1.0_c1.0

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list1[@]}"; do
    sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=B --meth=new --fix=$fix --a=1.0 --c=1.0 --logK=$logK |g" ../../submit_template.slm
    sbatch ../../submit_template.slm
  done
done
cd ../
### New KC-RPMD rates ###

cd ../
##############################
### System B a=1.0, c=1.0 ###
##############################

##############################
### System B a=0.1, c=0.1 ###
##############################
mkdir _sys_B_a0.1_c0.1
cd _sys_B_a0.1_c0.1

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list1[@]}"; do
    sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=B --meth=new --fix=$fix --a=0.1 --c=0.1 --logK=$logK |g" ../../submit_template.slm
    sbatch ../../submit_template.slm
  done
done
cd ../
### New KC-RPMD rates ###

cd ../
##############################
### System B a=0.1, c=0.1 ###
##############################

##############################
### System B a=0.1, c=5.0 ###
##############################
mkdir _sys_B_a0.1_c5.0
cd _sys_B_a0.1_c5.0

### New KC-RPMD rates ###
mkdir kcrpmd_new
cd kcrpmd_new
for fix in s y; do
  for logK in "${logK_list1[@]}"; do
    sed -i "s|python.*|python ../../1_kcrpmd_tst.py --sys=B --meth=new --fix=$fix --a=0.1 --c=5.0 --logK=$logK |g" ../../submit_template.slm
    sbatch ../../submit_template.slm
  done
done
cd ../
### New KC-RPMD rates ###

cd ../
##############################
### System B a=0.1, c=5.0 ###
##############################


