#!/usr/bin/env bash
ln -s ../nonUnitL14/interpoldata.dat blocks

for dsig in  1.1 1.2 1.3; do
    sed -e "28s/dd/$dsig/g" array.sh > array-$dsig.sh
    sed -e "28s/dd/$dsig/g" arrayNR.sh > arrayNR-$dsig.sh
    sbatch array-$dsig.sh
    sbatch arrayNR-$dsig.sh
done