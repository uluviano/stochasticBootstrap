#preamble
using DelimitedFiles
nRelevantScalars = 2
scalarGap = 1.3
include("utilities.jl");

tmpParam = readdlm("params.dat")
 stDim = tmpParam[1]|>Int


 DTfixed = tmpParam[2]|>Bool

#NR minimization by sector

nMins = 20

sectors = readdlm("sectorsJoined/sectors.dat", String)[:,1] |>unique
sectorsVec = Dict(sec =>[parse(Int, x ) for x in split(sec,'_')]  for sec in sectors)
lowTDataDump = []
for sec in sectors
    
    cd("sectorsJoined/lowT-$sec/")
   #jobID = readline("arrayJ")
    for i in 1:nMins
        try 
            nLinesStr = read(`wc -l putMin$i.dat`,String)
            if !occursin(r"^0 ",nLinesStr)
                dd =readdlm("putMin$i.dat")[end,:]
                inp = readdlm("putMin$i")
                sps = Int64.(inp[18:end,1] )
                push!(lowTDataDump,sectorSplitter(Float64.(dd), sps, 10^-8))
            end
        catch e
            print(e)
            @show sec, i
        end
    end
     cd("../../")
end


sectorList = [sectorToString(x[1]) for x in lowTDataDump]

@show sectors
@show unique(sectorList)

#NR minimization by ansatz

nMins = 20
nAnsatz = parse(Int, split(read(`wc -l ergodicNames`,String)," ")[1])
nOpsList = [1:nAnsatz;]

#lowTDataDump = []
for n in nOpsList
    
    cd("nOpsJoinedNR/nOps$n/")
   #jobID = readline("arrayJ")
    for i in 1:nMins
        try 
            nLinesStr = read(`wc -l putMin$i.dat`,String)
            if !occursin(r"^0 ",nLinesStr)
                dd =readdlm("putMin$i.dat")[end,:]
                inp = readdlm("putMin$i")
                sps = Int64.(inp[18:end,1] )
                push!(lowTDataDump,sectorSplitter(Float64.(dd), sps, 10^-8))
            end
        catch e
            print(" bad boys")
            @show n, i
        end
    end
     cd("../../")
end


@show nAnsatz

#clean up minima
uniqueSectors  = unique([x[1] for x in lowTDataDump])
dataBySector = Dict( sec => [] for sec in uniqueSectors)

filter!(
x->withinBounds(boundsFromSecnRelevantScalars(x[1],stDim;nRelevantScalars=nRelevantScalars, scalarGap=scalarGap, tol = 0.00),x[2]) , lowTDataDump
)
filter!(
x-> all(x[1] .!= 0) , lowTDataDump
)

for (sec, dd) in lowTDataDump
    push!(dataBySector[sec],dd )
end

relTol = 0.2
cleanMinima = vcat([[(sec,mm)  for mm in nonRedundantMinima(sec,minima, relTol )] for (sec, minima)  in dataBySector]...);
@show length(cleanMinima)

# generate minima checkBounds

        dirString = "checkingMinimaDependenceOnBoundaries"
        mkpath(dirString)
        cd(dirString)
  
		for (inputFileTag,(sec,dd)) in enumerate(cleanMinima)
			spinVector = spinVec(sec)
			 n = length(spinVector)
			 deltas = dd[2:n+1]
			 dsig = dd[1]
			



			dstar = maximum(spinVector) + stDim+1
			temp = 0.1

		#   generate input 


			movingOps = []
				operatorsInput = ["-1.	$dsig	$dsig	$dsig	0. \n"]

				operatorsGap = Dict(l=>0.0 for l in unique(spinVector))
			for (i,l) in enumerate(spinVector)

				bHi = dstar


						step =1

				if operatorsGap[l] == 0
					if l==0
						bLo =  (stDim-2) / 2 
						bHi = stDim
					else
						bLo =  l-2 + stDim
					end
					operatorsGap[l] = 0.11 #this can be changed to an increment to get a different gap for each op
					if l==2 && DTfixed
						step = 0
					end
					#dd = l-2 + stDim
				else
					#bLo =  l<=0 ? stDim : l-2 + stDim  + operatorsGap[l]
					bLo =  l<=0 ? 0 : l-2 + stDim  
					#dd = (bHi-bLo)rand()/2 + bLo
				end
				dd = deltas[i]
				if dd<bLo
					dd=bLo
				end

					
				push!(movingOps,step)
			push!(operatorsInput,"$l	$bLo	$bHi	$dd	$step \n")
			end
			
			gradientIndices = [1:n;][Bool.(movingOps)]
			gIstring = ""
			for x in gradientIndices
				gIstring*="$x "
			end

			inputFileStr = 
		join(["""blocks
			200
			0.42
			$dstar
			50000000
			1
			0
			$temp
			0.0001
			0.0001
			0.
			10000.
			1
			211
			$n
			$(n-1)
			$gIstring
			""", operatorsInput...] 
			)
			### PENDING ###
			# add name tag
			write("nbis$inputFileTag",inputFileStr)
			

		end

arrayLauncher=raw"""
#!/usr/bin/env bash
#
#SBATCH --job-name=bounds
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=uluviano@sissa.it
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-core=1
#
#SBATCH --mem-per-cpu=650mb
#
#SBATCH --array=1-"""*
"$(length(cleanMinima))"*
raw"""

#SBATCH --partition=regular2
#SBATCH --time=2:00:00
#SBATCH --output=%x.o%A-%a
#SBATCH --error=%x.e%A-%a
#

## YOUR CODE GOES HERE (load the modules and do the calculations)
## Sample code:

# Make sure it's the same module used at compile time
module load intel/2021.2


ln -s ../../nonUnitL14/interpoldata.dat blocks

# Run calculation
../lapackNewton< nbis${SLURM_ARRAY_TASK_ID} > nbis${SLURM_ARRAY_TASK_ID}.dat
		"""

			write("array-checkBounds.sh",arrayLauncher)
cd("../")
