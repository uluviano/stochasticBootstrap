#!/usr/bin/env bash
#
#SBATCH --job-name=jjj
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=uluviano@sissa.it
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-core=1
#
#SBATCH --mem-per-cpu=650mb
#
#SBATCH --array=1-220
#SBATCH --partition=regular2
#SBATCH --time=2:00:00
#SBATCH --output=%x.o%A-%a
#SBATCH --error=%x.e%A-%a
#

## YOUR CODE GOES HERE (load the modules and do the calculations)
## Sample code:

# Make sure it's the same module used at compile time
module load intel/2021


# Run calculation
./lapackNewton< lowT${SLURM_ARRAY_TASK_ID} > lowT${SLURM_ARRAY_TASK_ID}.dat

