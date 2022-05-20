#include("D:\\Github\\Projects\\Ongoing\\Relative_boolean\\bmodel_julia\\bmodel.jl")
# include("/mnt/d/Github/Projects/Ongoing/Relative_boolean/bm_julia_package/dependencyinstaller.jl")
include("/mnt/d/Github/Projects/Ongoing/Relative_boolean/bm_julia_package/bmodel.jl")
using Base.Threads

fileList = readdir()
topoFiles = String[]
for i in fileList
	if endswith(i, "topo")
		push!(topoFiles, i)
	end
end

#println(Threads.nthreads())

# topoFiles = ["EMT_RACIPE2.topo"]

Threads.@threads for topoFile in topoFiles
	y1 = @elapsed x = bmodel_reps(topoFile)
	println(topoFile, "-", y1)
	# y2 = @elapsed x = bmodel_reps(topoFile,100000, 1000, "Async", 0)
	# println(topoFile, " - ", y1, " and " , y2, " seconds.")
end