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
   dir=rootdir*"/Application"
   dirdata=dir*"/Data"
   dirresultsprelim=dir*"/Results_prelim"
   ## ######################### Functions #########################################
   include(dir*"/functions_decomp_app.jl")
end
## ######################### Loading Data #########################################
baseoption=1
nmarkets=34
alldata=Matrix(mydata(baseoption,nmarkets))
T=Int(maximum(alldata[:,2])); dY=Int(maximum(alldata[:,4]));
marketcol=15;
nmarkets=length(unique(alldata[:,marketcol]))

## Passing parameters for optimization to different procs
#Set default=3 to estimate the direct market shares (see Table 1).
@everywhere begin taun=0.01; default=0; alldata=$(alldata); marketcol=$(marketcol); end
## Optimization
println("Starting optimization")
@time begin Outcome=pmap(market->decomposition_market_app(market, alldata, marketcol, default), 1:nmarkets); end
## ######################### Collecting the results #########################################
println("Saving the results")
dA=length(Outcome[1][end])
ResultsDecomp=DataFrame(market_ids=Float64[],trip_ids=Int[], set_ids=Int[], product_ids=Int[], shares=Float64[], set_pr=Float64[])
for  t in 1:T, market in 1:nmarkets, prod in 1:dY, set in 1:dA
   push!(ResultsDecomp,[market t set prod Outcome[market][t][prod,set] Outcome[market][4][set]])
end
## Saving the output
CSV.write(dirresultsprelim*"/decomp_all_$(nmarkets)_$(default)_$(baseoption).csv", ResultsDecomp)
