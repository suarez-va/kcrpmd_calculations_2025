#!/bin/bash

cd _sys_A1_oldgamma

if [ -d "kcrpmd_ori" ]; then
  cd kcrpmd_ori
  for dir2 in _fix_*/; do
    cd "$dir2/libra_data" || continue
    pwd
    rm -r _itraj_*
    cd ../../
  done
  cd ../
fi

if [ -d "kcrpmd_new" ]; then
  cd kcrpmd_new
  for dir2 in _fix_*/; do
    cd "$dir2/libra_data" || continue
    pwd
    rm -r _itraj_*
    cd ../../
  done
  cd ../
fi
