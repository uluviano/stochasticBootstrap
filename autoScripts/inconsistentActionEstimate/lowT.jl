using DelimitedFiles
include("utilities.jl")

#get run global Params
tmpParam = readdlm("params.dat")
 stDim = tmpParam[1]|>Int


 DTfixed = tmpParam[2]|>Bool


  
sectors =  readdlm("sectors.dat")
sectorsVec = Dict(sec =>[parse(Int, x ) for x in split(sec,'_')]  for sec in sectors)



for sec in sectors
	secVec = sectorsVec[sec]
	ltpath = "lowT-$sec"
	mkpath(ltpath)
	cd(ltpath)
	n=sum(secVec)
	spinTuple = (-1,spinVec(secVec)...)
	dstar = maximum(spinTuple) + stDim+1
	hiBounds = tuple([dstar for _ in  1:length(spinTuple)]...)
	loBounds = tuple([l<=0 ? (stDim-2) ÷ 2 : l-2 + stDim for l in  [spinTuple...]]...)
	stepsVec = [1 for _ in 1:length(spinVec(secVec))]
	if DTfixed
		stepsVec[secVec[1]+1] = 0
	end
spinVector = spinTuple[2:end]
	steps = (0,stepsVec...)
				operatorsInputlowT = ["-1.	0	3	dd	0. \n"]

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
					# setting higher spin current
					if l==4
						step = 0
					end
					dd = l-2 + stDim
				else
					bLo =  l<=0 ? stDim : l-2 + stDim  + operatorsGap[l]
					dd = (bHi-bLo)rand()/2 + bLo
				end

				push!(movingOps,step)
			push!(operatorsInputlowT,"$l	$bLo	$bHi	dd	$step \n")
					
			end
			
			gradientIndices = [1:n;][Bool.(movingOps)]
			gIstring = ""
			for x in gradientIndices
				gIstring*="$x "
			end


			inputFileStr = 
		join([
    """blocks
200
0.42
$dstar
nSteps
1
0
0.tt
0.0001
0.0001
0.
10000.
10
211
nn 
$(n-2)
$gIstring
     """, operatorsInputlowT...] 
			)
			### PENDING ###
			# add name tag
	        write("inputFile",inputFileStr)
	
	launchFile = join( [raw"""
	#!/bin/bash
	#export LC_ALL="en_US.UTF-8"

	#This improved version uses variable size steps and allows for a very small step for the "fixed" operators in the exploratory phase

	ln -s ../../nonUnitL14/interpoldata.dat blocks

	putMinRank=1
	""" ,
	"nn=$n" ,
	raw"""

	while read p; do
		putMin=( $p)

	""",
	raw"""
	boundariesHi=($(awk '(NR>17){print $3}' inputFile))
boundariesLo=($(awk '(NR>17){print $2}' inputFile)) 
       steps=($(awk '(NR>17){print $5}' inputFile)) 
nvars=$(($nn - 1))
		""",
		raw"""

			echo ${putMin[@]}

            fileName=putMin${putMinRank}
            sed -e "8s/tt/00001/" -e "9s/MCstep/0.00001/" -e "13s/dens/1/" -e "5s/nSteps/500000000/" -e "15s/nn/$nn/" -e "16s/nvars/$nvars/" inputFile > $fileName
            row=$((18))
            for ii in ${!putMin[@]}; do
                row=$((18 + $ii))
                #echo $row
                step=1
                #correct boundaries
                #correct boundaries
                if (( $(echo "${putMin[${ii}]} < ${boundariesLo[${ii}]}" | bc -l)  )) ;then
                    echo ${putMin[${ii}]}
                    putMin[${ii}]=${boundariesLo[${ii}]}
                fi
                if (( $(echo "${putMin[${ii}]} > ${boundariesHi[${ii}]}" | bc -l)  )) ;then
                    echo ${putMin[${ii}]}
                    putMin[${ii}]=${boundariesHi[${ii}]}
                fi
                sed  -e "${row}s/bLo/${boundariesLo[${ii}]}/"  -e "${row}s/bHi/${boundariesHi[${ii}]}/" -e "${row}s/step/${steps[${ii}]}/" -e "${row}s/dd/${putMin[${ii}]}/" -i $fileName
            done

			#run

			echo "Runnning low T in all $nn directions"



		putMinRank=$(($putMinRank + 1))
	done <minima.dat
	"""]
	)
	write("test.sh",launchFile)

	cd("../")
end