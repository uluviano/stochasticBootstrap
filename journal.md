# 2022b07
Need to copy `testNR.sh` and `findNR.sh` to a folder and generate there the input files for the NR minimization of `2022a04-Ising-freetmuinu`

launched free Tmunu Nr minimizations by nnOps

launched free tmunu dsig = 2/7

# 2022b08
## TO DO
- [ ] add NR generation code to `sectors.jl`
- [ ] calculate spectrum for AT theories (code function)
  - [ ] decompose in SL(2,R) blocks
  - [ ] search for local minima nearby
- [ ] find the fate of 3d data, 
  - [ ] are there any stable minima?
  - [ ] if so, why NR and lowTMC give such different results?
- [ ] is there a real basin for the 2d ising?
- [ ] 

## Done:
- Analyzed free Tmunu 2d Ising data (see `/home/puercoerrante/fortranmc/dataAnalysisNotebooks/2022a04-2dIsing.ipynb`)
- Launched NR in dsig = 2/7 free Tmunu
- extended barriercheck to "" fixed Tmunu

# 2022b09
## TO DO
- [ ] add NR generation code to `sectors.jl`
- [ ] find the fate of 3d data, 
  - [ ] are there any stable minima?
  - [ ] if so, why NR and lowTMC give such different results?
- [ ] is there a real basin for the 2d ising?
- [ ] calculate spectrum for AT theories (code function)
  - [ ] decompose in SL(2,R) blocks
  - [ ] search for local minima nearby

## Done:

- found a bug in dealing with runaways (overflow lead to gammaN not going back after convergence).
- launched corrected barrier checks in dsig = 2/7 (free and fixed tmunu)
- analyzed them and realized that with the new convergence conditions we might be losing too much stuff.
- launched search for local minima in 3d Ising

# 2022b10
## TO DO
- [x] add NR generation code to `sectors.jl`
- [ ] find the fate of 3d data, 
  - [ ] are there any stable minima?
  - [ ] if so, why NR and lowTMC give such different results?
  - [ ] analyze 3d Ising MC data and launch nOps NR.
- [ ] is there a real basin for the 2d ising?
- [ ] calculate spectrum for AT theories (code function)
  - [ ] decompose in SL(2,R) blocks
  - [ ] search for local minima nearby

  ## done
- Corrected script with NR minimizations by sector and by fixed spectrum.
  - Solved issue with stDim consistency.
- Launched runs in 2d with other dsigma
  - 2022b10-dsig0.2312-dtFixedfalse
  - 2022b10-dsig0.2312-dtFixedtrue
- Launched 3d runs 
  - 2022b10-dsig0.5181489-dtFixedtrue
  - 2022b10-dsig0.5181489-dtFixedfalse

# 2022b11
## TO DO
- [ ] Close analysis of 2d data
  - [x] barrierCheck for dsig = 2/7, DTFixedtrue 
  - [x] barrierCheck for dsig = 2/7, DTFixedfalse 
    - [x] **Improvements**: add script that creates `array-barrierCheck.sh` in Julia
  - [ ] sectorsNR for dsig = 1/8, DTFixedfalse 
    - [ ] it seems like I already had the stuff
  - [x] sectorsNR for dsig = 0.2312, DTFixedfalse 
  - [x] sectorsNR for dsig = 0.2312, DTFixedtrue 
  - [x] nOpsNR for dsig = 0.2312, DTFixedfalse 
  - [x] nOpsNR for dsig = 0.2312, DTFixedtrue 
    - [ ] **Improvements**: add scripts for nOpsNR in autoScripts folder
``` 
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp findNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp arrayNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp testNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp minimaByNops-NR.sh ../2022b10-dsig0.2312-dtFixedfalse/ 
```


- [ ] find the fate of 3d data
  - [x] sectorsNR for dsig = ising, DTFixedfalse 
  - [x] sectorsNR for dsig = ising, DTFixedtrue 
  - [ ] nOpsNR for dsig = ising, DTFixedfalse 
  - [ ] nOpsNR for dsig = ising, DTFixedtrue 
  - [ ] are there any stable minima?
  - [ ] if so, why NR and lowTMC give such different results?
  - [ ] analyze 3d Ising MC data and launch nOps NR.
- [ ] is there a real basin for the 2d ising?
- [ ] calculate spectrum for AT theories (code function)
  - [ ] decompose in SL(2,R) blocks
  - [ ] search for local minima nearby

  ## done

# 2022b12
## TO DO
- [ ] Close analysis of 2d data
  - [ ] Analyze results and plot full picture (3 values of dsig and free vs fixed tmunu)
  - [ ] Need to find way to reduce number of minima
  - [ ] Understand bug in `ergodic.jl` that creates duplicate list in ergodicNames
    - [ ] **Improvements**: add scripts for nOpsNR in autoScripts folder
``` 
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp findNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp arrayNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp testNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp minimaByNops-NR.sh ../2022b10-dsig0.2312-dtFixedfalse/ 
```


- [ ] find the fate of 3d data
  - [ ] nOpsNR for dsig = ising, DTFixedfalse 
  - [ ] nOpsNR for dsig = ising, DTFixedtrue 
  - [ ] are there any stable minima?
  - [ ] if so, why NR and lowTMC give such different results?
  - [ ] analyze 3d Ising MC data and launch nOps NR.
- [ ] is there a real basin for the 2d ising?
- [ ] calculate spectrum for AT theories (code function)
  - [ ] decompose in SL(2,R) blocks
  - [ ] search for local minima nearby

  ## done

# 2022b14
## TO DO
- [ ] Close analysis of 2d data
  - [ ] Analyze results and plot full picture (3 values of dsig and free vs fixed tmunu)
  - [ ] Need to find way to reduce number of minima
  - [ ] Understand bug in `ergodic.jl` that creates duplicate list in ergodicNames
    - [ ] **Improvements**: add scripts for nOpsNR in autoScripts folder
``` 
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp findNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp arrayNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp testNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp minimaByNops-NR.sh ../2022b10-dsig0.2312-dtFixedfalse/ 
```


- [ ] find the fate of 3d data
  - [x] nOpsNR for dsig = ising, DTFixedfalse 
  - [x] nOpsNR for dsig = ising, DTFixedtrue 
  - [x] are there any stable minima?
  - [ ] if so, why NR and lowTMC give such different results?
  - [ ] analyze 3d Ising MC data and launch nOps NR.
- [ ] is there a real basin for the 2d ising?
- [ ] calculate spectrum for AT theories (code function)
  - [ ] decompose in SL(2,R) blocks
  - [ ] search for local minima nearby

  ## done
  - [x] nOpsNR for dsig = ising, DTFixedfalse 
  - [x] nOpsNR for dsig = ising, DTFixedtrue 
  - [ ] Started 3d analysis file

# 2022b15
## TO DO
- [x] Close analysis of 2d data
  - [ ] Analyze results and plot full picture (3 values of dsig and free vs fixed tmunu)
  - [ ] Need to find way to reduce number of minima
  - [ ] Understand bug in `ergodic.jl` that creates duplicate list in ergodicNames
    - [ ] **Improvements**: add scripts for nOpsNR in autoScripts folder
``` 
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp findNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp arrayNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp testNR.sh ../2022b10-dsig0.2312-dtFixedfalse/
[uluviano@login1 2022b10-dsig0.2312-dtFixedtrue]$ cp minimaByNops-NR.sh ../2022b10-dsig0.2312-dtFixedfalse/ 
```


- [ ] find the fate of 3d data
  - [x] nOpsNR for dsig = ising, DTFixedfalse 
  - [x] nOpsNR for dsig = ising, DTFixedtrue 
  - [x] are there any stable minima?
  - [ ] if so, why NR and lowTMC give such different results?
  - [ ] analyze 3d Ising MC data and launch nOps NR.
- [ ] is there a real basin for the 2d ising?
- [ ] calculate spectrum for AT theories (code function)
  - [ ] decompose in SL(2,R) blocks
  - [ ] search for local minima nearby

## done
  -Updated `stateOfTheArt.ipynb` with the 2d overview.
  -Analyzed 3d results around ising.

# 2022b16

## TODO:

## DONE:
  - analyzed 2d and 3d results together and did pdf.

# 2022b17

## TODO:
- [ ] find the fate of 3d data
  - [x] are there any stable minima?
    - [ ] if so, why NR and lowTMC give such different results?
- [x] 2d clarifications:
  - [x] Is 2d Ising stable along $\Delta_T$?

  ## DONE:
  - launched d=3 simulations
  - analyzed DT stability for dsig =1/8 in d=2

# 2022b18
- Low T analysis from 2d and 3d ising-like minima


# 2022b19
- took 2_1_1_1 minimum in d=3 (`/scratch/uluviano/mcboo/3d/2022b10-ising/checkingMinimaStabilityAlongdsig
`) and truncated to 2_1_1. 
	- NR at fixed dsig -> 2nd scalar decouples.
	- NR at fixed dsig (2nd scalar ~3.6 initial) -> 2nd scalar decouples. (AGAIN)
	- NR with free dsig (2nd scalar ~3.6 initial) -> 2nd scalar decouples. (AGAIN)

- adapted mathematicaTester to check lowT 2d evolutions
	- must analyze
# 2022b20 
- Discovered apparent stability was induced in many case due to wrong bounds. 
- Launched correct lowTMCs in d=2,3
  - Analyzed the results and found some minima stay around ising, some flow down to FT.

# 2022b21 

** See mail dsig = 2/7 **

# 2022b24
- v3 of quality criteria implementation gives acceptable results in dsig = 1/8 (local and non local) and in dsig = 2/7
- long MCs in d =2 alongg dsig are still non conclusive.

# 2022b25
**PENDING**
# 2022b26
- added function to plot spectrum nicely in `2022a04-2dIsing`
- added also cells to exclude bad endpoints and include "good" isolated points.

# 2022b28
- Polished 2d Plots. 
- Started final analysis of 3d data.

# 2022c01
- plotted minima with corresponding sector at each l_max: no clear correlation appears

# 2022c02
- added cell to remove redundant endpoints
- checked evolution of $\Delta_\epsilon$ in 2d Ising

# 2022c03
- agreed on the protocol + visualization exemplified in d=2 dsig

# 2022c10
- coded `getCleanMinima.jl` that goes into a directory with where `checkingMinimaDependenceOnBoundaries` has been computed and calculates relevant info:
  - `cleanMinima.dat`: all the stable, non redundant minima.
  - `specsToPlot.dat`: endpoints below a certain threshold and isolated minima within the lowest 0.05 lowest minima at each l_max.
  - `branch$ii.dat`: branches. Maybe not too useful.


# 2022c14
- Changed structure of `stateOfTheArt.ipynb`, now it's roughly: 
  - stDim
    - global analysis
    - interest cases
- started looking at DTvsDepsvsAct correlations.

# 2022c15
- caught up on missing runs in d=2,3. Launched final step of d=4.
- started actual sections in state of the art.

# 2022c16
- Checked the re-run of dsig =1/8 against the old results. They seem consistent enough.

# 2022c17
- solved graphical issues of global plot.


# 2022c18
- started chain computing function in state of the art
- created ergodic-oneOpPerSpin.jl that creates the input files for the specialized runs
- launched 
- - [ ] d=3 dsig = 1/2

  - [x] ergodic

  - [ ] NR

  - [ ] checkBounds

  - [ ] Analyze

- [ ] d=4 dsig = 1, 1.05, 1.15

  - [x] ergodic

  - [ ] NR

  - [ ] checkBounds

  - [ ] Analyze

- [ ] One op per spin re-do???

  - [x] ergodic

  - [ ] NR

  - [ ] checkBounds

  - [ ] Analyze

  # 2022c19
  - Launched NR for d=3, d=4 missing values, but found cluster issues.
  - Found a detail in module loading in nodes cn09-06, cn09-17,cn10-07,cn09-06, 
  - Started analyzing d=4 results, but found that in `/scratch/uluviano/PhD/mc-boo/4d/2022c12-dsig1.2-dtFixedtrue/checkingMinimaDependenceOnBoundaries` `nbis198.dat ` and`nbis199.dat ` have different dsig????

  # 2022c21
  - Corrected bug in `lowT.jl`.
  - added 3d scatter plots.
  - Launched Last bit of data generation.

  # 2022c22
  - Analysis of missing data points for 3d and 4d.
  - One op per spin in d=3,4.
  - Discovered possible bug in `boundaryCheck.jl` that rules out points close to the boundary. Analyzed `2022c18-oneOpPerSpin-dsig1.0-dtFixedtrue` and this seems to be the case.

  # 2022c23
  - Launched points in `3d/2022c23-inconsistent-dsig0.5*` with inconsistent setting of higher spin currents.
  - analyzed 3d ergodic runs to generate landscape.

  # 2022c24
  - relaunched MC part of grid search. Analyzed NR results, which look promising.
  - Studied the differences between the landscape for critical and multicritical Ising-like theories. 
    - How to choose "leading" scalar?

  # 2022c25
  - launched test for outliers in favorita
  - analyzed grid; might need higher density of points or re-running the start.
  - launched dsig @ ft test with tolerance for getting corner cases

  # 2022c26
  - launched test for outliers in favorita -> no strange features. added summary to tex
  - implemented branch explorer in stateOftheart

  # 2022c28
  - Launched initial part of lowTMC minimizations for showing convergence to FT
  - realized plotly sucks for latex :()

  # 2022c29
  - launched missing runs

  # 2022c30
  - retook FT analysis

  # 2022c31
  - finished FT analysis
  - relaunched last failed runs in d=3