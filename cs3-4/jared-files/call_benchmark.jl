using Distributed
using ClusterManagers
using BenchmarkTools
# using PyPlot
using DataFrames
using CSV

function run_the_benchmark(parray)

    tarray1 = zeros(length(parray))
    tarray2 = zeros(length(parray))

    for i in 1:length(parray)
        # addprocs(parray[i])
        addprocs(SlurmManager(parray[i]))
        include("benchmark_distributed_ieacs4.jl")
        # include("benchmark_distributed_eo2.jl")
        # println(mean(t1)[1])
        tarray1[i] = deepcopy(t1)
        tarray2[i] = deepcopy(t2)
        rmprocs(workers())
    end

    return tarray1, tarray2

end

parray = [1, 10, 50, 100, 200, 400, 800, 1440]

tarray1, tarray2 = run_the_benchmark(parray)

csvdata = DataFrame(np = parray, aep_time_sec = tarray1, jac_time_sec = tarray2)
CSV.write("par_scaling.csv", csvdata)

# plot(parray, tarray1/tarray1[1], label="AEP")
# legend()
# plt.show()