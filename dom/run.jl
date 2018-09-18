
@everywhere include("../RandomWalk.jl")

numtraj = 1000000
time = vcat(0, [2^(i/3) for i=-15:50])
Nvals = [2, 4, 6, 10, 12, 15, 20, 24, 30, 40, 50, 60, 80, 100]

for N in Nvals
    println("Begin N = $N")
    path = "/data/niels/sim/rw/dom/N_$N"
    mkpath(path)
    traj = RandomWalk.dom_random_walk(numtraj, 2*N, 70000, time)
    RandomWalk.savedata("$path/traj.h5", time, traj, 2*N)
end
