#This file generates the results for the 2nd set of simulations used in construction of Tables 3-10 in the online appendix
using Pkg
Pkg.activate(".")
using Distributed
addprocs(16) #Set number of available procs
@everywhere begin
   using Pkg
   Pkg.activate(".")
end
using CSV, DataFrames
@everywhere begin
   using LinearAlgebra
   using Random
   using Distributions, Statistics, StatsBase
   using Combinatorics
   using JuMP, KNITRO
## ######################### Dir ################################################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Replication_UCS",tempdir1)[end]]
dir=rootdir*"/Simulations/Recovering_FM/"
dirresults=dir*"Results_sim2/"
## ######################### Functions #########################################
include(dir*"sim2_functions.jl")
end
########################### Parameters #########################################
sampsize=2000 #Set the sample size sampsize either to 2000, 5000, 10000, 50000
model="core" #Set model either to "nested" or "core"
println(model)
println(sampsize)
## Distribution over the sets and choices
Pd,Pyd=dgpparam(model)
## Passing parameters for optimization to different procs
@everywhere begin
  Pd=$Pd; Pyd=$Pyd; sampsize=$sampsize
end
## Optimization
println("Starting optimization")
numMC=1000
@time Outcome=pmap(oneMC,1:numMC)

## Saving the results
M=DataFrame([Outcome[s][4] for s in 1:numMC],:auto)
F=vcat(DataFrame([Outcome[s][1][:] for s in 1:numMC],:auto),DataFrame([Outcome[s][2][:] for s in 1:numMC],:auto),DataFrame([Outcome[s][3][:] for s in 1:numMC],:auto))

CSV.write(dirresults*"M_$(model)_$(sampsize).csv", M)
CSV.write(dirresults*"F_$(model)_$(sampsize).csv", F)
