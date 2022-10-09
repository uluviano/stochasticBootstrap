#!/usr/bin/env bash
#
#SBATCH --job-name=boundaryCheck
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

# Run calculation
~/julia-1.6.5/bin/julia  boundaryCheck.jl

#need to test
cd checkingMinimaDependenceOnBoundaries
sbatch array-checkBounds.sh