#data pre-processing
function everyn(fileName::String,n::Int)
#`awk '!(NR % $nstr)' $fileName`
cc=`awk -v everyn=$n '!(NR % everyn )' $fileName`
outName="ev$(n)-$(fileName)"
open(outName, "w") do file                      
        write(file, read(cc))              
    end 
    return outName
end

function filterCols(fileName::String,inp::Vector{Int64})

cc=`awk -v x=$(inp[1]) -v y=$(inp[2]) -v z=$(inp[3]) '{print $x " " $y " " $z}' $fileName`
outName="rwready-$(fileName)"
open(outName, "w") do file                      
        write(file, read(cc))              
    end 
    return outName
end


function partitionVector(v::Vector,n::Int64)
    ll = length(v)
    println(ll)
    imax = Int(floor(ll/n))
    println(imax)
    return [v[1+(i-1)*(n):i*n] for i in 1:imax]
end

function sectorSplitter(frame::Vector{Float64}, spins::Vector{Int64}, tol::Float64)
   #do stuff heheee
   nOps = length(spins)
   ext = frame[1]
   act = frame[nOps+1]
   uniqueSpins = unique(spins[2:end])

   a=Dict(ll => Vector{Float64}[] for ll in uniqueSpins)
   for ii in 2:nOps
      if frame[ii + nOps] > tol && ( frame[ii] < 2*maximum(frame[deleteat!([1:nOps;],ii)]))
         push!(a[spins[ii]], frame[ii .+ [0, nOps]])
      end
   end

   for ll in uniqueSpins
      sort!(a[ll], by = x -> -x[2])
   end

   sector = [length(a[ll]) for ll in uniqueSpins]
   deltas = hcat([v[1] for ll in uniqueSpins for v in a[ll]]...)
   spins = hcat([v[2] for ll in uniqueSpins for v in a[ll]]...)

   return sector, hcat([ext], deltas,[act],spins)|>vec


end



function dataGrouper(frames::Matrix{Float64}, spins::Vector{Int64}, tol::Float64)
   #do stuff heheee
   nOps = length(spins)
   uniqueSpins = unique(spins[2:end])
   #maxSector = [count(==(ll), spins) for ll in uniqueSpins]
   fullClass = Dict{Vector{Int64}, Matrix{Float64}}()
   for ii in 1:size(frames)[1]
      ss, dd = sectorSplitter(frames[ii,:], spins, tol)
      if haskey(fullClass, ss)
         fullClass[ss] = vcat(fullClass[ss], dd)
      else
         fullClass[ss] =  dd
      end
   end

   return fullClass

end



function sectorToString(sector::Vector{Int64})
    result = ""
    for num in sector
        result *= "$(num)_"
    end
    return result[1:end-1]
end


#Contour plotting utilities 



function binActions(testData::Matrix{Float64},xs::AbstractRange{Float64},ys::AbstractRange{Float64},smooth::Bool=true,fillValue=Nothing,threshold =100)
 @show   nOps = size(testData)[2] ÷ 2
    binnedActions = [testData[(xx .< testData[:,1] .<= xx+xstep) .* (yy .< testData[:,2] .< yy + ystep),nOps + 1] for xx in xs for yy in ys ]
    function filledMin(list::Vector{Float64},fillValue)  
        if length(list)<threshold
            return fillValue
        else
            return min(list...)
        end
    end
    if fillValue == Nothing 
        fillValue = maximum(testData[:,nOps + 1])
    end
    binnedMinima = [filledMin(actions, fillValue) for actions in binnedActions]
    

    lx= length(ys)
    ly  = length(xs)
    matrixMinima = reshape(binnedMinima, (lx,ly) )
    if smooth
        for ii in 2:lx-1
            for jj in 2:ly-1 
                matrixMinima[ii,jj] = matrixMinima[ii,jj]/2 + ( matrixMinima[ii-1,jj] +matrixMinima[ii,jj+1] +matrixMinima[ii,jj-1] +matrixMinima[ii+1,jj] )/8
            end
        end
    end
    return matrixMinima
end

function nonRedundantMinima(meansFromCutoff::Vector, stdsFromCutoff::Vector,sigmaTol)
    trueMin = [1]
    for i in 2:length(meansFromCutoff)
        any(all(abs.(meansFromCutoff[j]-meansFromCutoff[i])/sqrt.(stdsFromCutoff[i].^2 .+ stdsFromCutoff[j].^2).<=sigmaTol) for j in trueMin) || push!(trueMin,i)
    end
    return trueMin
end
function nonzeroOPEs(meansFromCutoff::Vector,OPETol)
    nOps = size(meansFromCutoff[1])[2] ÷ 2
    trueMin = filter(ii -> all(abs.(meansFromCutoff[ii][nOps+2:end]) .> OPETol) , 1:length(meansFromCutoff))
    return trueMin
end
function converged(data::Vector,window)
    nOps = size(data[1])[2] ÷ 2
    trueMin = filter(ii -> minimum(data[ii][1:end - window,nOps+1]) <= minimum(data[ii][end - window+1:end,nOps+1]) , 1:length(data))
    return trueMin
end
function noCollision(meansFromCutoff::Vector,spins::Vector{Int64},tolerance)
    trueMin = Int64[]
    for ii in 1:length(meansFromCutoff)
        head = 2 
        good = true
        for spin in spins
            if spin>1
                deltas = sort(meansFromCutoff[ii][head:head+spin-1])
                a = deltas[1:end-1]
                b = deltas[2:end]
                if any( abs.(a .- b) .< tolerance)
                    good = false
                end
            end
            head += spin
        end
        good && push!(trueMin, ii)

    end
    return trueMin
end

function midpoint(x,y)
    (x+y)/2
end

struct Cell2d
    xlims::Tuple{Float64,Float64}
    ylims::Tuple{Float64,Float64}
    points::Matrix{Float64}
    cutoff::Int64

    function Cell2d(xlims,ylims,points,cutoff)
        if sum((xlims[1] .< points[:,1] .< xlims[2]) .* (ylims[1] .< points[:,2] .< ylims[2]))>cutoff

            #@show  ((xlims[1] .< points[:,1] .< xlims[2]) .* (ylims[1] .< points[:,2] .< ylims[2]))
            new(xlims,ylims, points[((xlims[1] .< points[:,1] .< xlims[2]) .* (ylims[1] .< points[:,2] .< ylims[2])),:], cutoff)
        else
            new(xlims,ylims,[ NaN NaN NaN ],cutoff)
        end
    end

 end

function quarter(cell::Cell2d)
    midpointx = midpoint(cell.xlims...)
    midpointy = midpoint(cell.ylims...)
    xintervals = [(cell.xlims[1], midpointx ), ( midpointx, cell.xlims[2] )]
    yintervals = [(cell.ylims[1], midpointy ), ( midpointy, cell.ylims[2] )]
    [Cell2d(xx, yy , cell.points, cell.cutoff) for xx in xintervals for yy in yintervals]
end


spinVec(spins) = [ 2(i-1) for i in 1:length(spins) for j in 1:spins[i]]
ref3dIsingFull=vcat([[0.5181489  -1  1], [1.412625 0 1.1063962062036898], [3 2 0.1063658384978176], [5.022665 4 0.004771493775999999], [3.82968 0 0.002810272144], [6.42065 4 3.822807039999999e-6], [7.028488 6 0.00024779797056000006], [5.50915 2 0.00011182005025000001], [7.38568 4 5.6522685025000005e-6], [6.8956 0 5.384624399999999e-7], [7.2535 0 2.6244000000000003e-8], [7.0758 2 2.2781529e-7]]...)
ref3dIsing=vcat([[0.5181489  -1  1], [1.412625 0 1.1063962062036898], [3 2 0.1063658384978176], [5.022665 4 0.004771493775999999], [3.82968 0 0.002810272144], [6.42065 4 3.822807039999999e-6], [7.028488 6 0.00024779797056000006], [5.50915 2 0.00011182005025000001], [6.8956 0 5.384624399999999e-7]]...)
ref2dIsing = [
1/8 -1  1
1	0	1/4
2	2	1/64
4	4	9/40960
4	0	1/4096
5	4	1/65536
6	6	25/3670016
6	2	9/2621440
7	6	1/1310720
8	8	15527/57579405312
8	4	25/234881024
8	0	81/1677721600
9	8	1125/30064771072
9	0	1/1073741824
10	10	251145/20882130993152
10	6	15527/3685081939968
10	2	45/30064771072
11	10	227/120259084288
11	2	1/21474836480
12	12	18598401/32317945275219968
12	8	251145/1336456383561728
12	4	46581/786150813859840
12	0	625/13469017440256
13	12	454541/4644337115725824
13	4	1125/492581209243648
13	0	1/429496729600
14	14	2638578917/91484645086776524800
14	10	18598401/2068348497614077952
14	6	452061/171066417095901184
14	2	388175/211317338765524992
15	14	208917/40250921669623808
15	6	227/1970324836974592
15	2	225/1970324836974592
16	16	11593243155/7762080058562616229888
16	12	2638578917/5855017285553697587200
16	8	167385609/1323743038473009889280
16	4	6278625/76637754858963730432
16	0	241087729/3315387916083573817344
17	16	25880735367/92325954088916305838080
17	8	454541/76092819304051900416
]
ref3dFT=vcat([[0.5  -1  1], [1 0 2.], [3 2 0.09375], [5 4 0.00427246], [7 6 0.000220299], [9 8 0.0000119861]]...)


function smartBin(testData::Matrix{Float64},xinitial,yinitial, threshold::Int64,cutoff::Int64)
    #xlims = extrema(testData[:,1])
    #ylims = extrema(testData[:,2])
    #bigCell = Cell2d(xlims,ylims,testData)

    xboundaries = [(xinitial[i],xinitial[i+1]) for i in 1:length(xinitial)-1]
    yboundaries = [(yinitial[i],yinitial[i+1]) for i in 1:length(yinitial)-1]
    #reduce
    cellList = [Cell2d(xlim,ylim,testData,cutoff) for xlim in xboundaries for ylim in yboundaries]
    i =1
    while any(length(cell.points) >1000 for cell in cellList) && i<threshold
        @show i+=1
        cellList = vcat([length(cell.points) >threshold ? quarter(cell) : [cell] for cell in cellList]...)
    end

    #create grid
    xs = unique([cell.xlims[1] for cell in cellList])|>sort!
    ys = unique([cell.ylims[1] for cell in cellList])|>sort!
    xstep = minimum(xs[2:end] - xs[1:end-1])
    ystep = minimum(ys[2:end] - ys[1:end-1])
    xvec = range(xs[[1,end]]...,step = xstep)
    yvec = range(ys[[1,end]]...,step = ystep)
    resultMatrix = zeros((length(xvec),length(yvec)))
    for ix in 1:length(xvec), iy in 1:length(yvec)
    cell = filter(i -> (i.xlims[1] <= xvec[ix] <= i.xlims[2] ) && (i.ylims[1] <= yvec[iy] <= i.ylims[2] ), cellList )[1]
    resultMatrix[ix,iy] = minimum(cell.points[:,3])[1]
    end
    return xvec, yvec, resultMatrix
end

#plotting

function grayScaleCustom(nS, grayCutoff, minVal,maxVal)
    step = (maxVal - minVal)/nS
    [[ i/nS, "rgb($(grayCutoff + i*(255-grayCutoff) ÷ nS), $(grayCutoff + i*(255-grayCutoff)÷nS),$(grayCutoff + i*(255-grayCutoff)÷nS) )"] for i in 0:nS]


end

function connectedPairs(xs,ys,eps,subsDim)
    distMatrix = [sum( ((x[1:subsDim]-y[1:subsDim]) ./ y[1:subsDim]).^2 ) |>sqrt for x in xs, y in ys ] 
    findall(x->x<eps, distMatrix)
end

colorsMinima = [
"black",
"blueviolet",
"brown",
"burlywood",
"cadetblue",
"chartreuse",
"chocolate",
"coral",
"cornflowerblue",
"cornsilk",
"crimson",
"cyan",
"darkblue",
"darkcyan",
"darkgoldenrod",
"darkgray",
"darkgrey",
"darkgreen",
"darkkhaki",
"darkmagenta",
"darkolivegreen",
"darkorange",
"darkorchid",
"darkred",
"darksalmon",
"darkseagreen",
"darkslateblue",
"darkslategray",
"darkslategrey",
"darkturquoise",
"darkviolet",
"deeppink",
"deepskyblue",
"dimgray",
"dimgrey",
"dodgerblue",
"firebrick",
"floralwhite",
"forestgreen",
"fuchsia",
"gainsboro",
"ghostwhite",
"gold",
"goldenrod",
"gray",
"grey",
"green",
"greenyellow",
"honeydew",
"hotpink",
"indianred",
"indigo",
"ivory",
"khaki",
"lavender",
"lavenderblush",
"lawngreen",
"lemonchiffon",
"lightblue",
"lightcoral",
"lightcyan",
"lightgoldenrodyellow",
"lightgray",
"lightgrey",
"lightgreen",
"lightpink",
"lightsalmon",
"lightseagreen",
"lightskyblue",
"lightslategray",
"lightslategrey",
"lightsteelblue",
"lightyellow",
"lime",
"limegreen",
"linen",
"magenta",
"maroon",
"mediumaquamarine",
"mediumblue",
"mediumorchid",
"mediumpurple",
"mediumseagreen",
"mediumslateblue",
"mediumspringgreen",
"mediumturquoise",
"mediumvioletred",
"midnightblue",
"mintcream",
"mistyrose",
"moccasin",
"navajowhite",
"navy",
"oldlace",
"olive",
"olivedrab",
"orange",
"orangered",
"orchid",
"palegoldenrod",
"palegreen",
"paleturquoise",
"palevioletred",
"papayawhip",
"peachpuff",
"peru",
"pink",
"plum",
"powderblue",
"purple",
"red",
"rosybrown",
"royalblue",
"saddlebrown",
"salmon",
"sandybrown",
"seagreen",
"seashell",
"sienna",
"silver",
"skyblue",
"slateblue",
"slategray",
"slategrey",
"snow",
"springgreen",
"steelblue",
"tan",
"teal",
"thistle",
"tomato",
"turquoise",
"violet",
"wheat",
"white",
"whitesmoke",
"yellow",
"yellowgreen"
]


## Functions for minima connection

function compatibleSectors(va::Vector, vb::Vector, maxDif::Int64; onlyContained = true, excludeIdentical = true )
    #handle identical sectors exception:
    if excludeIdentical
        va==vb && return false
    end
    #ensure hierarchy
    length(va)>length(vb) && return false
    v1  = copy(va)
    v2  = copy(vb)


    minDim,maxDim = length(v1),length(v2)
    #minDim>maxDim && return false
    for _ in 1:(maxDim-minDim)
        push!(v1,0)
    end
    if onlyContained
        return sum( abs.(v1- v2)) <=maxDif &&  all(v1 .<= v2 )
    else
        return sum( abs.(v1- v2)) <=maxDif #&&  all(v1 .<= v2 )
    end
end


function vecToSpinDict(frame::Vector{Float64}, spinMults::Vector{Int64})
   #do stuff heheee
   spins = spinVec(spinMults)
   nOps = length(spins)

   ext = frame[1]
   act = frame[nOps+2]
   uniqueSpins = unique(spins)

   spinDict = Dict(ll => Float64[] for ll in uniqueSpins)
   for ii in 1:nOps
      push!(spinDict[spins[ii]], frame[ii + 1])

   end

   return spinDict


end

function maxDistBetweenMinima(va::Vector{Float64}, vb::Vector{Float64}, seca::Vector{Int64}, secb::Vector{Int64})
    maxl = 2*min(length(seca), length(secb)) -2 

    da = vecToSpinDict(va, seca)
    db = vecToSpinDict(vb, secb)
distances = [[ [ abs((x - y)/sqrt(x^2 + y^2)) for x in da[ll] , y in db[ll]][i,i] for i in 1:min(length(da[ll]),length(db[ll])) ] 
        for ll in 0:2:maxl
        ]
    maximum(
       vcat(distances...        )
    )
    

end


         function sectorToString(sector::Vector{Int64})
            result = ""
            for num in sector
         result *= "$(num)_"
            end
            return result[1:end-1]
         end



function areMinimaCompatible(va::Vector{Float64}, vb::Vector{Float64}, seca::Vector{Int64}, secb::Vector{Int64}, relTol::Float64)
   maxDistBetweenMinima(va,vb,seca,secb) < relTol

end

function nonRedundantMinima(sec::Vector{Int},minima::Vector{Any},relTol::Float64)  
    if length(minima)>1
        n = sum(sec)
        distMatrix = [maxDistBetweenMinima(m1,m2,sec,sec) for m1 in minima, m2 in minima]  
        minimaToKeep=[1:length(minima);]
        for i in 1:length(minima)
            closeCalls = findall(x->x<relTol, distMatrix[i,i+1:end])
            for j in closeCalls
                if minima[i][n+2] <= minima[j][n+2]
                    filter!(x->x!=j,minimaToKeep)
                else
                    filter!(x->x!=i,minimaToKeep)
                end
            end


        end
        unique!(minimaToKeep)
        return minima[minimaToKeep]
    else
        return minima
    end
end



struct Chain
    list::Vector
end


function containerVector(vecList::Vector)
    targetDim = maximum(length.(vecList))
    tmp=deepcopy(vecList)
    for v in tmp
        neededZeros = targetDim - length(v)
        for _ in 1:neededZeros
            push!(v, 0)
        end
    end
    maximum(hcat(tmp...), dims=2)[:,1]
end


function weakestLink(chain::Vector{Tuple{Vector{Int64}, Vector{Float64}}})
    clen = length(chain)
    [
maxDistBetweenMinima(
chain[i][2]  ,
chain[i+1][2]  ,
chain[i][1]  ,
chain[i+1][1]
)
for i in 1:(clen-1)
]
end



function withinBounds(bounds,dd)
    all(bounds[1] .<= dd[1:length(bounds[1])]) && all(bounds[2] .>= dd[1:length(bounds[2])])
end

function boundsFromSec(sec,stDim; omitUnitarity=true, gap = 1, tol = 0.01)
    spinVector = spinVec(sec)
			dstar = maximum(spinVector) + stDim+1
				operatorsGap = Dict(l=>0 for l in unique(spinVector))
            lower= [-1.0]
            upper= [dstar for i in 0:length(spinVector)]
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
					if omitUnitarity
						bLo -= 10 * tol
					end
					operatorsGap[l] = gap #this can be changed to an increment to get a different gap for each op
				else
					bLo =  l<=0 ? stDim : l-2 + stDim  + operatorsGap[l]
					dd = (bHi-bLo)rand()/2 + bLo
				end
                
					
			push!(lower, bLo)
			end
            return(lower .+ tol ,upper)
end

function recursiveNodeConnector(chain::Vector{Tuple{Vector{Int64}, Vector{Float64}}},
goodPairs::Vector,
 roots::Vector{Tuple{Vector{Int64}, Vector{Float64}}})
 
 if chain[end] in roots || length(chain) > 6
    return Chain(chain)
 else
    return [ recursiveNodeConnector(vcat(chain,[(s2,d2)]),goodPairs, roots) for ((s1,d1),(s2,d2)) in goodPairs if (s1,d1) == chain[end] ]
 end
end

function cleanNestedLists(chainList)
    while typeof(chainList) != Vector{Chain}
        chainList = vcat(chainList...)
    end
    return [chain.list for chain in chainList]
end


function atSpectrum(g,dStar)

    dmn(m,n) = (1+g)n^2/(1-g)+ (1-g)m^2/(1+g)
    lmn(m,n) = 2m*n
    #return dd, spins
end


function boundsFromSecScalarGap(sec,stDim;  gap = 1, tol = 0.01)
    spinVector = spinVec(sec)
			dstar = maximum(spinVector) + stDim+1
				operatorsGap = Dict(l=>0 for l in unique(spinVector))
            lower= [-1.0]
            upper= [dstar for i in 0:length(spinVector)]
			for l in spinVector

				bHi = dstar


						step =1

				if operatorsGap[l] == 0
					if l==0
						bLo =  (stDim-2) / 2 

						bHi = gap
					else
						bLo =  l-2 + stDim
					end
					operatorsGap[l] = 1 #this can be changed to an increment to get a different gap for each op
				else
					bLo =  l<=0 ? gap : l-2 + stDim  + operatorsGap[l]
					dd = (bHi-bLo)rand()/2 + bLo
				end
                
					
			push!(lower, bLo)
			end
            return(lower .+ tol ,upper)
end

function boundsFromSecnRelevantScalars(sec,stDim;  nRelevantScalars = 1, tol = 0.01)
    spinVector = spinVec(sec)
			dstar = maximum(spinVector) + stDim+1
				operatorsGap = Dict(l=>0 for l in unique(spinVector))
            lower= [-1.0]
            upper= [dstar for i in 0:length(spinVector)]
				scalarCount=0
			for l in spinVector

				bHi = dstar


						step =1

				if operatorsGap[l] == 0
					if l==0
						bLo =  (stDim-2) / 2 
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
						bLo = (stDim-2) / 2 
						bHi = stDim
						scalarCount+=1
					end
					dd = (bHi-bLo)rand()/2 + bLo
				end
                
					
			push!(lower, bLo)
			end
            return(lower .+ tol ,upper)
end

function boundsFromSecnRelevantScalars(sec,stDim;  nRelevantScalars = 1, scalarGap = 1.0, tol = 0.01)
    spinVector = spinVec(sec)
			dstar = maximum(spinVector) + stDim+1
				operatorsGap = Dict(l=>0 for l in unique(spinVector))
            lower= [-1.0]
            upper= [dstar for i in 0:length(spinVector)]
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
                
					
			push!(lower, bLo)
			end
            return(lower .+ tol ,upper)
end