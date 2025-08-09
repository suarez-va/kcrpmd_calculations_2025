#!/bin/bash

for dir in _sys_*/; do
  cd $dir
  pwd
  rm -r _itraj_*
  rm -r kappa_data/
  cd ../
done

