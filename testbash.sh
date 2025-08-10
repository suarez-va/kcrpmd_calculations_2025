#!/bin/bash

#python kcrpmd_thermalization.py --sys=1 --method=1 --fix=s --K0=0.00285 --leps=0.015 --hw=0
a_list=(0.1 0.1 0.1 0.3 1.0 1.0 1.0 1.0 1.0 1.0)
K0_list=(9.55e-06 9.55e-05 3.02e-04 5.37e-04 9.55e-04 1.70e-03 3.02e-03 4.47e-03 6.46e-03 9.33e-03)
leps_list=(-1.24e-02 -1.33e-02 -1.43e-02 -1.52e-02 -1.62e-02)


for i in "${!K0_list[@]}"; do
  K0=${K0_list[i]}
  if [[ $i -lt 3 ]]; then
    continue
  fi
done
