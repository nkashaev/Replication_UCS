#This file generates the results for the 1st set of simulations used in construction of Tables 1-2 in the online appendix
using Pkg
Pkg.activate(".")
using Distributed
addprocs(16)#Set number of available procs
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
dir=rootdir*"/Simulations/Recovering_Sets/"
dirresults=dir*"Results_sim1/"
## ######################### Functions #########################################
include(dir*"sim1_functions.jl")
end
## ######################### Parameters #########################################
sampsize=2000 #Set the sample size sampsize either to 2000, 5000, 10000, 50000
model="nested" #Set model either to "nested" or "core"
println(model)
println(sampsize)
## Distribution over the sets and choices
Pd,Pyd=dgpparam(model)
## Passing parameters for optimization to different procs
@everywhere begin
  tau=0.01; Pd=$Pd; Pyd=$Pyd; sampsize=$sampsize
end
## Optimization
println("Starting optimization")
numMC=1000
@time begin Outcome=pmap(oneMC,1:numMC); end
## Preparing the output of MC simulation
Ms1=DataFrame([Outcome[s][4] for s in 1:numMC],:auto)
Ms2=DataFrame([Outcome[s][8] for s in 1:numMC],:auto)
## Saving the results
CSV.write(dirresults*"Ms1_$(model)_$(sampsize).csv", Ms1)
CSV.write(dirresults*"Ms2_$(model)_$(sampsize).csv", Ms2)

