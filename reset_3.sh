#!/bin/bash

for dir in _sys_3*/; do
  cd "$dir/libra_data"
  pwd
  rm -r _itraj_*
  cd ../../
done

