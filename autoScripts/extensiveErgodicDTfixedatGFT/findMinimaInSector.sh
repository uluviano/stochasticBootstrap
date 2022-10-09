#!/usr/bin/env bash
#
#SBATCH --job-name=prev-secc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=uluviano@sissa.it
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-core=1
#
#SBATCH --mem-per-cpu=850mb
#
#SBATCH --partition=regular1
#SBATCH --time=6:00:00
#SBATCH --output=%x.o%A-%a
#SBATCH --error=%x.e%A-%a
#

## YOUR CODE GOES HERE (load the modules and do the calculations)
## Sample code:

# Make sure it's the same module used at compile time

module load intel/2021


ln -s ../../../nonUnitL14/interpoldata.dat blocks

	outname=secc
	echo $outname

	n=$(sed -e "1s/_/ /g" <<< $outname | awk 'BEGIN{tot = 0} {for(i=1;i<=NF; i++){ tot = tot+$i }} END{print tot}')
	echo "nOps:"
	echo $n

	sed -e "17s/nn/$n/" ../../find_min.f90 > f$n.f90
	ifort -O3 f$n.f90 -o f$n
	echo "finding minima"
	
	./f$n <  ../tmp-erg-$outname.dat > ranked-$outname.dat 
	sort -nrk1 ranked-$outname.dat | head -n20 | awk '{for(i=2; i<NF-1; i++) printf $i " "; print $(NF-1)}' > minima.dat;
	# launching lowT
	cd lowT-$outname/;  
	echo "running low T"
	. test.sh &>launch.log   
	sed -e "3s/jjj/$outname/" ../../array.sh> array.sh
	#cp ../../lowTMC .

	sbatch array.sh
