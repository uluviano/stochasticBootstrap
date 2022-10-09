#!/bin/bash
#export LC_ALL="en_US.UTF-8"


ln -s ../nonUnitL14/interpoldata.dat blocks
while read n; do 
echo $n
	sed -e"27s/jjj/$n/g" -e"3s/jjj/$n/" ergodic.sh > n$n.sh
	sbatch n$n.sh
 done < ergodicNames