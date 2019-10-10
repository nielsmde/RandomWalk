using Distributed 

@everywhere include("../RandomWalk.jl")


numtraj = 1000000
time = vcat(0, [2^(i/3) for i=-15:50])
# Nvals = [2, 4, 6, 10, 12, 15, 20, 24, 30, 40, 50, 60, 80, 100]
Nvals = [4, 10, 20]
biases = [-2.0, 2.0, -1.5, 1.5, -1.0, 1.0,-0.6, 0.6]

for N in Nvals, b in biases
    println("Begin N = $N, b = $b")
    path = "/data/niels/sim/rw/biased/N_$(N)_b_$b"
    mkpath(path)
    traj = RandomWalk.biased_random_walk(numtraj, 2*N, 300000, time, bias=b)
    RandomWalk.savedata("$path/traj.h5", time, traj, 2*N)
end
