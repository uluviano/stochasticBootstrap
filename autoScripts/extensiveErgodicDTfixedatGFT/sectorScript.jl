#using Pkg
#Pkg.add("DelimitedFiles")
using DelimitedFiles
  

function sectorSplitter(frame::Vector{Float64}, spins::Vector{Int64}, tol::Float64)
   #do stuff heheee
   nOps = length(spins)
   ext = frame[1]
   act = frame[nOps+1]
   uniqueSpins = unique(spins[2:end])

   a=Dict(ll => Vector{Float64}[] for ll in uniqueSpins)
   for ii in 2:nOps
      if frame[ii + nOps] > tol
         push!(a[spins[ii]], frame[ii .+ [0, nOps]])
      end
   end

   for ll in uniqueSpins
      sort!(a[ll], by = x -> -x[2])
   end

   sector = [length(a[ll]) for ll in uniqueSpins]
   deltas = hcat([v[1] for ll in uniqueSpins for v in a[ll]]...)
   spins = hcat([v[2] for ll in uniqueSpins for v in a[ll]]...)

   return sector, hcat([ext], deltas,[act],spins)


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


function splitter(listData::Vector{String}, step, spins::Vector{Int64}, fileName::String)
   for ff in listData
      println("importing data...")
      @time fullErgodicData = readdlm("$(ff)")
      println("done ")


      #add for here. Split it up and use write
      nFrames = length(fullErgodicData[:,1])

      nBatches = ceil(nFrames/step)|>Int
      for i in 1:nBatches
         @show 1 + step*(i-1):step*(i)
         if i== nBatches
            @time mm = dataGrouper(fullErgodicData[1 + step*(i-1):end, :], spins, 10^-8)
         else
            @time mm = dataGrouper(fullErgodicData[1 + step*(i-1):step*(i), :], spins, 10^-8)
         end

         #clearing
         #fullErgodicData = nothing

         function sectorToString(sector::Vector{Int64})
            result = ""
            for num in sector
         result *= "$(num)_"
            end
            return result[1:end-1]
         end

         for sector in keys(mm)
            sector_string = sectorToString(sector)
            open("$fileName-$sector_string.dat", "a") do io
         writedlm(io,  mm[sector])
            end

         end
   end
      end
end


dirName = "sectorsJoined"
mkpath(dirName)

fnames = readdlm("ergodicNames")
#fnames = ["n$n" for n in 4:22]
for ff in fnames
   fname = "$ff.dat"
   #input
   spins=readdlm("$ff")[15:end,1] |> Vector{Int}

	@show fname
	@show spins
	splitter([fname], 10000, spins, "$dirName/erg")
end
