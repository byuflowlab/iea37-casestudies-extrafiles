using ClusterManagers

addprocs(SlurmManager(4), t="00:01:00")

hosts = []
pids = []
for i in workers()
        host, pid = fetch(@spawnat i (gethostname(), getpid()))
        println(host)
        push!(hosts, host)
        push!(pids, pid)
end


# The Slurm resource allocation is released when all the workers have
# exited
for i in workers()
        rmprocs(i)
end

println(hosts, pids)