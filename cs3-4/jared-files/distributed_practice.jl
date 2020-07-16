using Distributed
using ClusterManagers


addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))
# addprocs(0)
# hosts = []
# pids = []
# for i in workers()
#         host, pid = fetch(@spawnat i (gethostname(), getpid()))
#         println(host)
#         push!(hosts, host)
#         push!(pids, pid)
# end
t1 = time()
c = @sync @distributed (+) for i = 1:10

    sleep(1)
    b = 1

end
t2 = time()

# The Slurm resource allocation is released when all the workers have
# exited
for i in workers()
        rmprocs(i)
end

println("ran in ", t2-t1)
println(c)
