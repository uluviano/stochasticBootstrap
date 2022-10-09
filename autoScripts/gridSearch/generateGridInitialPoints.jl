using DelimitedFiles
include("utilities.jl")

#Global Params
for stDim in 2:4

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
	DTfixed = true
	inputFileList = []
        dirString = "2022c23-gridInitialPoints-d$stDim"
        mkpath(dirString)
        cd(dirString)
        for ff in files
            cp("../$ff", ff)
        end
  
		inputFileTag = 1
		lmax = 6
		extraops = [
			0
			0
			2
			4

		]
		dsig0 = (stDim-2)/2 + 0.005
		dsig1 = (stDim-2)/2 + 0.2
		deps0 = (stDim-2) - 0.1
		deps1 = (stDim-2) + 0.9
		for dsig in dsig0:0.01:dsig1, deps in deps0:0.1:deps1
			@show dsig, deps
			
			#keep one in three

			#check good spin partition hierarchy: monotonically decreasing n(l)
			good=true
			nold = sum(extraops.==0)
			for lcheck in 2:2:lmax
				nnew = sum(extraops.==lcheck)
				if nold < nnew
					good=false
					break
				end
				nold = nnew
			end
			good || continue
			@show  inputFileTag

			spinVector =  vcat([0:2:lmax;],extraops)
			n = length(spinVector)
			
			dstar = maximum(spinVector) + stDim+1
			temp = 0.0001

		#   generate input 


			operatorsInput =     ["-1.	$((stDim-2) / 2)	$stDim	$dsig	0.1 \n"]
			operatorsInputlowT = ["-1.	$((stDim-2) / 2)	$stDim	$dsig	1 \n"]

			operatorsGap = Dict(l=>0 for l in unique(spinVector))
			movingOps=[]
			for l in spinVector

				bHi = dstar


				step =1

				dd = l-2 + stDim
				if operatorsGap[l] == 0
					if l==0
						bLo =  (stDim-2) / 2 
						bHi = stDim
						dd = deps
					else
						bLo =  l-2 + stDim
					end
					operatorsGap[l] = 1 #this can be changed to an increment to get a different gap for each op
					if l==2 && DTfixed
						step = 0
					end
				else
					#bLo =  l<=0 ? stDim : l-2 + stDim  + operatorsGap[l]
					# to see flow properly we remove the gaps
					bLo =  l-2 + stDim 
					dd = (bHi-bLo)rand()/2 + bLo
				end

				push!(movingOps,step)
			push!(operatorsInputlowT,"$l	$bLo	$bHi	$dd	$step \n")
					
			push!(operatorsInput,"$l	$bLo	$bHi	$dd	$step \n")
			end
			
			gradientIndices = [1:n;][Bool.(movingOps)]
			gIstring = "0 "
			for x in gradientIndices
				gIstring*="$x "
			end


			inputFileStr = 
		join(["""blocks
			200
			0.42
			$dstar
			200000000
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
			$(n)
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
			100000
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