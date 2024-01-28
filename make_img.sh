#!/bin/bash

#declare -a lstDir
lstDir=("A" "B" "C" "D" "E" "F")

#declare -a patchDir
patchDir=("X-Disp" "X-Loading" "XY-Rot" "Y-Disp" "Y-Loading")

for patch in "${lstDir[@]}"; do
  for testcase in "${patchDir[@]}"; do
    X:/Programme/Paraview/bin/pvpython.exe pvPatchTest.py "patch$patch-$testcase"
    #echo "patch$patch-$testcase"
  done
done

read