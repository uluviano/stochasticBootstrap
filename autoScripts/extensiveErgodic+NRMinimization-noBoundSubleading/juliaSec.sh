#!/usr/bin/env bash
#
#SBATCH --job-name=sectors
#SBATCH --mail-type=ALL
#SBATCH --mail-user=uluviano@sissa.it
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-core=1
#
#SBATCH --mem-per-cpu=1400mb
#
#SBATCH --partition=regular1
#SBATCH --time=1:00:00
#SBATCH --output=%x.o%A-%a
#SBATCH --error=%x.e%A-%a
#

## YOUR CODE GOES HERE (load the modules and do the calculations)
## Sample code:

# Make sure it's the same module used at compile time

module load intel/2021

# Run calculation
~/julia-1.6.5/bin/julia  sectorScript.jl

cp lowT.jl sectorsJoined/
cp lapackNewton sectorsJoined/
cp params.dat sectorsJoined/
cp utilities.jl sectorsJoined/

cd sectorsJoined

wc -l  erg-*.dat |head -n -1> namesTmp.dat


while read fname; do
	nLines=$(awk '{print $1}' <<< $fname)
	dname=$(awk '{print $2}' <<< $fname)
	nMin=120000
	if [[ "$fname" != *"_0"*  ]] && [ $nLines -ge $nMin ]; then
		echo "sector $dname: "
		denFact=$(($nLines/$nMin))
		echo $denFact
		awk -v nn=$denFact '!(NR%nn)' $dname > tmp-$dname
	 	sed -e "1s/-/ /g" -e "1s/\./ /g" <<< $dname| awk '{print $2}' >> sectors.dat
	fi

done < namesTmp.dat





#run Julia script for generating input files
~/julia-1.6.5/bin/julia  lowT.jl


# cd into each directory in sectors.dat
while read fname; do
	cd lowT-$fname
	sed -e "3s/secc/$fname/" -e "29s/secc/$fname/" ../../findMinimaInSector.sh > ff.sh
	sbatch ff.sh
	cd ../

done < sectors.dat

# launch script that computes minima and launches array

# profit!
