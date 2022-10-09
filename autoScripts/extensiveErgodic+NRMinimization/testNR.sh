#!/bin/bash
#export LC_ALL="en_US.UTF-8"

#This improved version uses variable size steps and allows for a very small step for the "fixed" operators in the exploratory phase

ln -s ../../nonUnitL14/interpoldata.dat blocks

putMinRank=1
nn=nnn
while read p; do
	putMin=( $p)

	boundariesHi=($(awk '(NR>17){print $3}' inputFile))
boundariesLo=($(awk '(NR>17){print $2}' inputFile)) 
       steps=($(awk '(NR>17){print $5}' inputFile)) 
nvars=$(($nn - 1))

		echo ${putMin[@]}

		fileName=putMin${putMinRank}
		sed -e "8s/tt/00001/" -e "9s/MCstep/0.00001/" -e "13s/dens/10/" -e "5s/nSteps/500000000/" -e "15s/nn/$nn/" -e "16s/nvars/$nvars/" inputFile > $fileName
		row=$((18))
		for ii in ${!putMin[@]}; do
			row=$((18 + $ii))
			#echo $row
			step=1
			#correct boundaries
			#correct boundaries
			if (( $(echo "${putMin[${ii}]} < ${boundariesLo[${ii}]}" | bc -l)  )) ;then
				echo ${putMin[${ii}]}
				putMin[${ii}]=${boundariesLo[${ii}]}
			fi
			if (( $(echo "${putMin[${ii}]} > ${boundariesHi[${ii}]}" | bc -l)  )) ;then
				echo ${putMin[${ii}]}
				putMin[${ii}]=${boundariesHi[${ii}]}
			fi
			sed  -e "${row}s/bLo/${boundariesLo[${ii}]}/"  -e "${row}s/bHi/${boundariesHi[${ii}]}/" -e "${row}s/step/${steps[${ii}]}/" -e "${row}s/dd/${putMin[${ii}]}/" -i $fileName
		done


		#run

		echo "Runnning low T in all $nn directions"



	putMinRank=$(($putMinRank + 1))
done <minima.dat
