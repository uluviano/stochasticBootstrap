using DelimitedFiles
include("utilities.jl")

#Global Params
stDim = 3

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
			6 => 0.6
			8 => 0.7
			10 => 0.9
		]
	)
end

if stDim == 4                                                                                      
	tVec = Dict(
		[
			4 => 0.2
			6 => 0.6
			8 => 0.7
			10 => 0.9
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
 


#get filelist
files = readdir()
for dsig in ([1.0:0.1:1.3;]./2)
    for DTfixed in [ true, false]
inputFileList = []
        dirString = "2022c18-oneOpPerSpin-dsig$dsig-dtFixed$DTfixed"
        mkpath(dirString)
        cd(dirString)
        for ff in files
            cp("../$ff", ff)
        end
  
		inputFileTag = 1
		iCont = -1
		for lmax in [4,6,8,10 ]
			
			#keep one in three


			 spinVector =  [0:2:lmax;]
			 n = length(spinVector)
			
			dstar = maximum(spinVector) + stDim+1
			temp = tVec[lmax]

		#   generate input 


				operatorsInput = ["-1.	$dsig	$dsig	$dsig	0. \n"]
				operatorsInputlowT = ["-1.	$dsig	$dsig	$dsig	0. \n"]

				operatorsGap = Dict(l=>0 for l in unique(spinVector))
				movingOps=[]
			for l in spinVector

				bHi = dstar


						step =1

				if operatorsGap[l] == 0
					if l==0
						bLo =  (stDim-2) / 2 
						bHi = stDim
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
