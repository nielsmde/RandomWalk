# Simulation module to perform random walk simulations

This Julia package allows for performing random walk simulations of 
stochastic jumps on a periodic grid, consisting of alternating domains with
varying mobility.

## Usage

Refer to the official [Julia website](https://julialang.org/) on how to setup a Julia environment on your system.

An exemplary Julia script how to use the packge is give in the file 
`dom/run.jl`, which should be run like this, to perform simulation for
multiple domain sizes N in parallel:
    
    julia -p auto run.jl

The produced trajectories are stored as HDF5 files in the repsective subdirecories.

The top level function to perform simualitons is 

    RandomWalk.dom_wandom_walk(numtraj, N, time, rate_1=1e-2, rate_2=1)
    
The parameters `numtraj` and `N` define the number of *particles* to simulate
and the total size of the grid,respectively; accordingly, the size of each domain is `N/2`.
The time steps of the resulting trajectory are defined by the parameter `time`,
which can be spaced unevenly, if desired. 
The mobility in both domains can be specified by the two parameters `rate_1/2`.
