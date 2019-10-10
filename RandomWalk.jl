module RandomWalk

using Distributed

using Interpolations
using HDF5

export dom_random_walk, biased_random_walk

struct Trajectory
    positions::Array{Int, 2}
    times::Array{Float64, 1}
    box::Int
end

function translate(s::Array{Int, 1}, step::Array{Int, 1}, N::Int)
    p  = s + step
    for i = 1:3
        if p[i] > N
            p[i] = p[i] - N
        elseif p[i] < 1
            p[i] = N - p[i]
        end
    end
    return p
end

function next_step(traj::Trajectory, step::Int, rates)
    # choose random number and calc wait times
    z = rand(2 * size(traj.positions, 2))
    tw = - 1 ./ rates(traj.positions[step, :]) .* log.(1 .- z)    
    # select direction of shortest tw
    t = minimum(tw)  
    i = indexin([t], tw)[1]
    d = ceil.(Int, i/2)
    p1 = traj.positions[step, :]
    p1[d] += 2 * (i % 2) - 1
    # append new position to trajectory
    traj.times[step + 1] = traj.times[step] + t
    traj.positions[step + 1, :] = p1
end

function random_walk(steps::Int, timemax::Float64, rates, N; dim::Int=3)
    traj = Trajectory(Array{Int}(undef, steps, dim), Array{Float64}(undef, steps), N)
    traj.times[1] = 0.0
    traj.positions[1, :] = ceil.(Int, rand(dim) * N)
    for i = 1:steps-1
        if traj.times[i] > timemax
            # when maximum time is reached, fill array and end loop
            traj.times[i + 1:end] .= traj.times[i] + 1.0
            # println("Reached maxtime in $i steps.")
            break            
        end
        next_step(traj, i, rates)
    end
    if traj.times[end] < timemax
        println("Maxtime was not reached: times[end] = $(traj.times[end])")
    end
    return traj
end


function to_uniform_time(time::AbstractArray{Float64, 1}, tr::Trajectory)
    pos = Array{Int}(undef, length(time), size(tr.positions, 2))
    j = 1
    for (i, t) in enumerate(time)
        while (j < length(tr.times)) & (tr.times[j+1] <= t)
            j += 1
        end
        # pos[i, :] = tr.positions[tr.times .<= t, :][end, :]
        pos[i, :] = tr.positions[j, :]
    end
    return pos
end

function savedata(fname::String, time, traj, box)
    # traj has shape (times, 3, numtraj)
    h5open(fname, "w") do file
        write(file, "box", box)
        write(file, "time", time)
        g = HDF5.g_create(file, "positions")
        for i = 1:length(time)
            g["$(i-1)"] = traj[i, :, :]
        end
    end
end

function dom_random_walk(numtraj::Int, N::Int, steps::Int, time::AbstractArray{Float64,1}; rate_1=1e-2, rate_2=1.0)#::Array{Float64, 3}
    timemax = maximum(time)
   
    traj = @distributed (a, b) -> cat(a, b, dims=3) for i = 1:numtraj
        dom_rates(pos::Array{Int, 1}) = Bool(sum(mod1.(pos, N) .> N/2) % 2) ? rate_1 : rate_2
        tr = random_walk(steps, timemax, dom_rates, N)
        to_uniform_time(time, tr)
    end
    return traj
end

function biased_random_walk(numtraj::Int, N::Int, steps::Int, time::AbstractArray{Float64, 1}; rate_1=1e-2, rate_2=1.0, bias=0.25)
    timemax = maximum(time)

    traj = @distributed (a, b) -> cat(a, b, dims=3) for i = 1:numtraj
        
        function biased_rates(pos::Array{Int, 1})
            x = mod1.(pos, N)[1]
            if (x == 1) | (x == N)
                b = exp(bias)
            elseif (x == (N / 2)) | (x == (N / 2 + 1))
                b = exp(-bias)
            else
                b = 1
            end
            # b = 1 + bias * cos(2 * pi * (pos[1] - 0.5) / N)
            k = Bool(sum(x .> N/2) % 2) ? rate_1 : rate_2
            [1/b, b] .* k
        end

        tr = random_walk(steps, timemax, biased_rates, N, dim=1)
        to_uniform_time(time, tr)
    end
    return traj

end

end
