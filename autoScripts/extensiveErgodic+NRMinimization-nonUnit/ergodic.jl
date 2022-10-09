using DelimitedFiles
include("utilities.jl")

#Global Params
stDim = 2

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


		  function powerset(x::Vector{T}) where T
			result = Vector{T}[[]]
			for elem in x, j in eachindex(result)
				push!(result, [result[j] ; elem])
			end
			result
		end
 


#get filelist
files = readdir()
for dsig in [-0.4]
    for DTfixed in [ true, false]
inputFileList = []
        dirString = "2022d15-dsig$dsig-dtFixed$DTfixed"
        mkpath(dirString)
        cd(dirString)
        for ff in files
            cp("../$ff", ff)
        end
  
		inputFileTag = 1
		iCont = -1
		for lmax in [4,6,8,10 ], extraops in unique(powerset([ll for ll in 0:2:(lmax-2) for _ in 1:min(4,(lmax-ll)÷2) ]))
			
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
			@show iCont, inputFileTag

			iCont +=1
			iCont%3!=0 && continue 

			 spinVector =  vcat([0:2:lmax;],extraops)
			 n = length(spinVector)
			
			dstar = maximum(spinVector) + stDim+1
			temp = tVec[lmax]

		#   generate input 


				operatorsInput = ["-1.	$dsig	$dsig	$dsig	0. \n"]
				operatorsInputlowT = ["-1.	$dsig	$dsig	$dsig	0. \n"]
				unitTol = 0.9

				operatorsGap = Dict(l=>0 for l in unique(spinVector))
				movingOps=[]
			for l in spinVector

				bHi = dstar


						step =1

				if operatorsGap[l] == 0
					dd = l-2 + stDim
					if l==0
						bLo =  (stDim-2) / 2  - unitTol
						bHi = stDim
						dd = dsig
					else
						bLo =  l-2 + stDim- unitTol
					end
					operatorsGap[l] = 1 #this can be changed to an increment to get a different gap for each op
					if l==2 && DTfixed
						step = 0
					end
				else
					bLo =  l<=0 ? stDim : l-2 + stDim  + operatorsGap[l]- unitTol
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
			0
			1
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
			0
			1
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
