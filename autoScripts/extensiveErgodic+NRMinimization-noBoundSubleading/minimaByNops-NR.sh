#!/usr/bin/env bash

mkdir nOpsJoinedNR

cp lapackNewton nOpsJoinedNR/

cd nOpsJoinedNR



while read name; do 
	n=$(sed -e"1s/n//" <<<$name)
	echo $n
	mkdir nOps$n
	cd nOps$n
	sed -e "3s/secc/$n/" -e "29s/secc/$n/" ../../findNR.sh > ff.sh
	sbatch ff.sh
	cd ../

 done < ../ergodicNames
