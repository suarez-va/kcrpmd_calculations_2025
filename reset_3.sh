#!/bin/bash

for dir in _sys_*/; do
  cd "$dir/libra_data"
  pwd
  rm -r _itraj_*
  cd ../../
done

