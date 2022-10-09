using DelimitedFiles

#Global Params
for stDim in 4:4

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
	DTfixed = false
	inputFileList = []
        dirString = "2022e06-GFTsearch-d$stDim"
        mkpath(dirString)
        cd(dirString)
        for ff in files
            cp("../$ff", ff)
        end
  
			

		for dsig in [1.1, 1.2, 1.3]
		inputFileTag = 1

			for lmax in 4:2:10
				fullSpectrum = [( ll+2dsig+n, ll) for n in 2:2:lmax+2 for ll in  0:2:lmax-n +2]
				#check good spin partition hierarchy: monotonically decreasing n(l)
				for lengthExtra in 0:length(fullSpectrum)
					spec = fullSpectrum[1:lengthExtra]
					good=true
					extraops = [x[2] for x in spec]
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
					deltasVector = vcat([0:2:lmax;],[x[1] for x in spec])
					n = length(spinVector)
					
					dstar = maximum(spinVector) + stDim+1
					temp = 0.0001

				#   generate input 


					operatorsInput =     ["-1.	$((stDim-2) / 2)	$stDim	$dsig	0 \n"]
					operatorsInputlowT = ["-1.	$((stDim-2) / 2)	$stDim	$dsig	0 \n"]

					operatorsGap = Dict(l=>0 for l in unique(spinVector))
					movingOps=[]
					for (opi,l) in enumerate(spinVector)

						bHi = dstar


						step =1

						dd = l + 2dsig
						if operatorsGap[l] == 0
							if l==0
								bLo =  (stDim-2) / 2 
								bHi = stDim
								dd = 2dsig
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
							dd = deltasVector[opi]
						end

						push!(movingOps,step)
					push!(operatorsInputlowT,"$l	$bLo	$bHi	$dd	$step \n")
							
					push!(operatorsInput,"$l	$bLo	$bHi	$dd	$step \n")
					end
					
					gradientIndices = [1:n;][Bool.(movingOps)]
					gIstring = " "
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
					write("lowT$inputFileTag-dsig$dsig",inputFileStr)
					

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
					nameInput = "n$inputFileTag-dsig$dsig"
					write(nameInput,inputFileStr)

					push!(inputFileList, nameInput)
					inputFileTag += 1
				end

			end
		end
        writedlm("ergodicNames",inputFileList)
        writedlm("params.dat",[stDim, DTfixed])

		cd("../")
end