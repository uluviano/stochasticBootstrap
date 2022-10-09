using DelimitedFiles
include("utilities.jl")

function cleanSpectrum!(spec::Vector{Int})
    x0 = spec[1]
    for i in 2:length(spec)
        if spec[i] > x0
            spec[i] = x0
        end
        x0 = spec[i]
    end
    return spec
end

#Global Params
stDim = 3
nRelevantScalars = 2
scalarGap = 1.3

if stDim == 2
	tVec = Dict(
		[
			4 => 0.3
			6 => 0.7
			8 => 0.8
			10 => 0.9
		]
	)
end                                                                                                
                                                                                                   
if stDim == 3                                                                                      
	tVec = Dict(
		[
			4 => 0.2
			6 => 0.5
			8 => 0.6
			10 => 0.8
		] 
	)
end


		  function powerset(x::Vector{T}) where T
			result = Vector{T}[[]]
			for elem in x, j in eachindex(result)
				push!(result, [result[j] ; elem])
			end
			result
		end
 

dsigsToPlot = [
     0.503
 0.505
 0.508
 0.51
 0.513
 0.514
 0.515
 0.516
 0.517
 0.518
 #0.5185
 0.519
 #0.5195
 0.52
 0.521
 0.522
 0.523
 0.524
 0.527
 0.53
]

#get filelist
files = readdir()
#for dsig in [0.534]
#for dsig in vcat([0.5195,0.5185],0.513:0.001:0.521)
for dsig in dsigsToPlot
    for DTfixed in [ true, false]
inputFileList = []
        dirString = "2022g06-dsig$dsig-dtFixed$DTfixed-ONexperiments"
        mkpath(dirString)
        cd(dirString)
        for ff in files
            cp("../$ff", ff)
        end

		#improved spectrum ansatz

		#obtained from previous runs
baseSpecs=[[2, 1, 1], [2, 1, 1], [2, 2, 1], [2, 1, 1, 1], [2, 2, 1, 1], [3, 2, 1, 1], [2, 3, 1, 1], [3, 3, 1, 1], [3, 2, 1, 1, 1], [2, 1, 1, 1, 1], [2, 2, 1, 1, 1], [2, 2, 2, 1, 1], [5, 3, 1, 1, 1], [3, 2, 1, 2, 1, 1], [4, 3, 1, 2, 1, 1], [4, 3, 1, 1, 1, 1], [3, 1, 2, 1, 1, 1], [3, 2, 2, 2, 1, 1], [3, 3, 2, 2, 1, 1], [2, 2, 2, 1, 1, 1], [5, 2, 1, 1, 1, 1], [3, 4, 1, 2, 1, 1], [2, 2, 2, 2, 1, 1], [3, 1, 1, 1, 1, 1], [4, 2, 2, 1, 1, 1], [3, 3, 1, 1, 1, 1], [3, 2, 2, 3, 1, 1], [5, 3, 1, 1, 1, 1], [4, 2, 1, 2, 1, 1], [4, 2, 1, 1, 1, 1], [4, 1, 2, 1, 1, 1], [3, 3, 2, 1, 1, 1], [3, 2, 3, 2, 1, 1], [3, 2, 1, 1, 1, 1], [3, 3, 1, 2, 1, 1], [3, 1, 2, 2, 1, 1], [4, 2, 2, 3, 1, 1], [4, 3, 2, 1, 1, 1], [4, 2, 2, 2, 1, 1], [4, 3, 2, 2, 2, 1], [3, 2, 2, 1, 1, 1], [5, 2, 2, 1, 1, 1], [3, 2, 3, 1, 1, 1], [4, 3, 2, 2, 1, 1], [2, 2, 2, 2, 2, 1], [4, 2, 1, 3, 1, 1], [2, 2, 1, 1, 1, 1], [3, 3, 1, 1, 2, 1], [2, 2, 1, 1, 2, 1], [4, 2, 3, 1, 1, 1], [4, 4, 2, 2, 1, 1], [4, 3, 1, 2, 2, 1], [5, 2, 2, 2, 1, 1], [3, 3, 2, 2, 2, 1], [4, 3, 3, 2, 1, 1], [3, 1, 3, 1, 1, 1], [3, 3, 2, 1, 2, 1], [2, 2, 1, 2, 2, 1], [2, 3, 2, 2, 1, 1], [4, 2, 1, 2, 2, 1], [3, 4, 2, 2, 1, 1], [4, 4, 1, 1, 1, 1], [3, 3, 3, 2, 1, 1], [2, 2, 1, 2, 1, 1], [3, 4, 1, 1, 1, 1], [3, 4, 3, 1, 1, 1], [4, 2, 2, 2, 2, 1], [3, 3, 3, 1, 1, 1], [3, 3, 1, 3, 1, 1], [3, 2, 1, 2, 2, 1], [3, 3, 2, 3, 1, 1], [4, 3, 3, 1, 1, 1], [3, 2, 1, 3, 1, 1], [3, 3, 3, 1, 2, 1], [3, 2, 4, 2, 1, 1], [5, 2, 1, 2, 1, 1], [3, 2, 2, 2, 2, 1], [4, 4, 2, 1, 1, 1], [2, 4, 2, 2, 1, 1], [2, 1, 2, 2, 1, 1], [4, 2, 3, 2, 1, 1], [3, 4, 2, 3, 1, 1], [5, 2, 1, 2, 2, 1]]
#clean spectra (enforce hierarchy)

cleanSpectrum!.(baseSpecs)

		inputFileTag = 1
		iCont = -1
		for spec in baseSpecs
			
			@show spec, inputFileTag


			 spinVector =  spinVec(spec)
			 n = length(spinVector)
			 lmax = maximum(spinVector)
			
			dstar = maximum(spinVector) + stDim+1
			temp = tVec[lmax]
			@show lmax dstar temp

		#   generate input 


				operatorsInput = ["-1.	$dsig	$dsig	$dsig	0. \n"]
				operatorsInputlowT = ["-1.	$dsig	$dsig	$dsig	0. \n"]

				operatorsGap = Dict(l=>0 for l in unique(spinVector))
				movingOps=[]
				scalarCount=0
			for l in spinVector

				bHi = dstar


						step =1

				if operatorsGap[l] == 0
					if l==0
						bLo =  0.9 
						bHi = stDim
						scalarCount+=1
					else
						bLo =  l-2 + stDim
					end
					operatorsGap[l] = 1 #this can be changed to an increment to get a different gap for each op
					if l==2 && DTfixed
						step = 0
					end
					dd = l-2 + stDim
				else
					bLo =  l<=0 ? stDim : l-2 + stDim  + operatorsGap[l]
					if l==0 && scalarCount < nRelevantScalars
						bLo = scalarGap 
						bHi = stDim
						scalarCount+=1
					end
					dd = (bHi-bLo)rand()/2 + bLo
				end
				 

				push!(movingOps,step)
			push!(operatorsInputlowT,"$l	$bLo	$bHi	dd	$step \n")
					
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
			10
			211
			$n
			$(n-1)
			$gIstring
			""", operatorsInputlowT...] 
			)
			### PENDING ###
			# add name tag
			write("lowT$inputFileTag",inputFileStr)
			

			inputFileStr = 
		join(["""blocks
			200
			0.42
			$dstar
			200000000
			1
			0
			$temp
			0.001
			0.
			10000.
			1000
			211
			$n
			""", operatorsInput...] 
			)
			### PENDING ###
			# add name tag
			nameInput = "n$inputFileTag"
			write(nameInput,inputFileStr)

            push!(inputFileList, nameInput)
			inputFileTag += 1

		end
        writedlm("ergodicNames",inputFileList)
        writedlm("params.dat",[stDim, DTfixed])

		cd("../")
	end 
end
