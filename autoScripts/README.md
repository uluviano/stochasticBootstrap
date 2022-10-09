# README #

This folder contains Julia scripts in order to generate comprehensive tests to be run in a cluster with SLURM as a scheduler. Each subfolder contains slightly different test cases.  *TODO: refactor as a single routine*

# To generate the scripts #
In the appropriate folder, edit `ergodic.jl` in line 5 and set `stDim` to the spacetime dimension of interest and in line 42 `dsig` to the values of  $\Delta_\sigma$ to be studied. The run in a terminal
`julia ergodic.jl` and the folders to be copied to the cluster will be generated.

# To run in the cluster #
Copy to the cluster and run the following scripts in order (i.e. wait for the previous step to finish):
1) `launchErgodic.sh`
2) `juliaSec.sh` and `findNR.sh`
3) `boundaryCheck.sh`

Then the data can be downloaded and analyzed as shown in `../resultsForPaper.ipynb`