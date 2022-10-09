# README #

This is v0.1 of the MonteCarloBootstrap project.

### Contains  ###

* Fortran Source code.
* Examples


### How do I get set up? ###

To run this package you need:

* An up-to-date fortran compiler (ifort if you have an intel processor, otherwise gfortran works too)
* (Optional:) Julia binary. 

### Usage ###

In order to explore the landscape you need:

* An interpolator data file (`interpoldata.dat` is provided for examples).
* An `input` file with the spectrum characteristics and MC params. (See `Implementation_summary.md` for details.)

After compiling `gfortran -O3 -llapack MCLAPACK.f90 -o a` you just need to type `./a < input > output.dat` and the evolution of the MC will be stored in `output.dat`.
(`newton.f90` is similar but uses Newton-Rhapson to find local minima deterministically.)

In the folder `autoScripts` some comprehensive tests can be generated.

