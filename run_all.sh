#!/bin/bash
nstep=1000
dt=1.0
a_list=(0.002 0.05 0.1 0.1 0.1 0.1 0.1)
K0_list=(0.0000095 0.0001710 0.0005320 0.0009500 0.0016911 0.0030021 0.0094909)

for ((i=0; i<${#a_list[@]}; i++))
do
    sed -i "s/python.*/kcrpmd_thermalization.py --nsteps $nstep --dt $dt --a ${a_list[$i]} --temp 300.0 --K0 ${K0_list[$i]} --bq 3.0 --eps 0.0 --systype B --fix y /g" submit_template.slm
    echo "python kcrpmd_thermalization.py --nsteps $nstep --dt $dt --a ${a_list[$i]} --temp 300.0 --K0 ${K0_list[$i]} --bq 3.0 --eps 0.0 --systype B --fix y /g"
    
done

