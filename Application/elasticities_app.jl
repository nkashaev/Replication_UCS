#Combining the results of decomposition with attributes, prices, and instruments
using Pkg
Pkg.activate(".")
using CSV, DataFrames
using Statistics, LinearAlgebra
using JuMP, KNITRO
########################### Dir ################################################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Replication_UCS",tempdir1)[end]]
dir=rootdir*"/Application"
dirdata=dir*"/Data"
dirresultsprelim=dir*"/Results_prelim"
dirtabfig=dir*"/Tables_and_figures"
######################## Parameter used in reading the data
nmarkets=34 
baseoption=1
T=3; dY=5; marketcol=15; default=0;
########################### Functions #########################################
include(dir*"/functions_elasticities_app.jl")
################################################################################
# Estimated beta
data_naiv=data_gen(31,baseoption,nmarkets, 3, 7);
beta_naiv, se_naiv=gmmbeta(data_naiv,4,1)


data=data_gen(31,baseoption,nmarkets, default, 7);
beta_decom31, se_decom31=gmmbeta(data,4,0)

data=data_gen(30,baseoption,nmarkets, default, 7);
beta_decom30, se_decom30=gmmbeta(data,4,0)

Table7=DataFrame(betase=[],naiv=[],set31=[],set30=[])
push!(Table7,["beta" beta_naiv beta_decom31 beta_decom30])
push!(Table7,["se" se_naiv se_decom31 se_decom30])
CSV.write(dirtabfig*"/Table13oa.csv",Table7)

############################
## Elasticities for Market 1
Elas=zeros(dY,4)
#Average elsasticity
data=data_gen(31,baseoption,nmarkets, 3, 7);
dataFC=data_gen(31,baseoption,nmarkets, default, 7);
dataNoQ=data_gen(30,baseoption,nmarkets, default, 7);
m=1;
for i in 1:dY
    indt=(data.product_ids.==i).*(data.trip_ids.==1).*(data.market_ids.==m)
    indtFC=(dataFC.product_ids.==i).*(dataFC.trip_ids.==1).*(dataFC.market_ids.==m)
    indtNoQ=(dataNoQ.product_ids.==i).*(dataNoQ.trip_ids.==1).*(dataNoQ.market_ids.==m)
    Elas[i,1]=beta_naiv*data.price[indt][1]*(1.0 - data.shares[indt][1])*(data.shares[indt][1].!=0)
    Elas[i,3]=beta_decom31*dataFC.price[indtFC][1]*(1.0 - dataFC.shares[indtFC][1])*(dataFC.shares[indtFC][1].!=0)
    Elas[i,4]=beta_decom30*dataNoQ.price[indtNoQ][1]*(1.0 - dataNoQ.shares[indtNoQ][1])*(dataNoQ.shares[indtNoQ][1].!=0)
    elFC=beta_decom31*dataFC.price[indtFC][1]*(1.0 - dataFC.shares[indtFC][1])*dataFC.shares[indtFC][1]/data.shares[indt][1]
    elNoQ=(beta_decom30*dataNoQ.price[indtNoQ][1]*(1.0 - dataNoQ.shares[indtNoQ][1])*dataNoQ.shares[indtNoQ][1]/data.shares[indt][1])*(dataNoQ.shares[indtNoQ][1].!=0)
    Elas[i,2]=elFC*dataFC.set_pr[indtFC][1]+elNoQ*dataNoQ.set_pr[indtNoQ][1]
end
CSV.write(dirtabfig*"/Table2.csv",DataFrame(Elas,:auto))

## Median Elasticities
Elas=zeros(dY,3)
#Naive
data=data_gen(31,baseoption,nmarkets, 3, 7);
Fcmarkets=unique(data.market_ids)
for i in 1:dY
    temp=zeros(length(Fcmarkets))
    for m in 1:length(Fcmarkets)
        indt=(data.product_ids.==i).*(data.trip_ids.==1).*(data.market_ids.==Fcmarkets[m])
        temp[m]=beta_naiv*data.price[indt][1]*(1.0 - data.shares[indt][1])
    end
    Elas[i,1]=median(temp)
end

#FC
data=data_gen(31,baseoption,nmarkets, default, 7);
Fcmarkets=unique(data.market_ids)
for i in 1:dY
    temp=zeros(length(Fcmarkets))
    for m in 1:length(Fcmarkets)
        indt=(data.product_ids.==i).*(data.trip_ids.==1).*(data.market_ids.==Fcmarkets[m])
        temp[m]=beta_decom31*data.price[indt][1]*(1.0 - data.shares[indt][1])
    end
    Elas[i,2]=median(temp)
end

#No Q
data=data_gen(30,baseoption,nmarkets, default, 7);
Fcmarkets=unique(data.market_ids)
for i in 1:dY-1
    temp=zeros(length(Fcmarkets))
    for m in 1:length(Fcmarkets)
        indt=(data.product_ids.==i).*(data.trip_ids.==1).*(data.market_ids.==Fcmarkets[m])
        temp[m]=beta_decom30*data.price[indt][1]*(1.0 - data.shares[indt][1])
    end
    Elas[i,3]=median(temp)
end


CSV.write(dirtabfig*"/Table14oa.csv",DataFrame(Elas,:auto))
