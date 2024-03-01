#Combining the results of decomposition with attributes, prices, and instruments
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Replication_UCS",tempdir1)[end]]
using Pkg
cd(rootdir*"/Application/Env 1.10")
Pkg.activate(".")
using CSV, DataFrames
using Statistics, LinearAlgebra
using Plots

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
include(dir*"/functions_dataandest_app.jl")
################################################################################
######## Summary Statistics
alldata=mydata(baseoption,nmarkets);
Sampsize=length(unique(alldata.household_code))
SumStat=zeros(3,5)
temp=[Vector(alldata[findfirst(alldata.ziphhn.==m),18:20]) for m in 1:nmarkets]
for k in 1:3
    SumStat[k,1]=mean([temp[m][k] for m in 1:nmarkets])
    SumStat[k,2]=median([temp[m][k] for m in 1:nmarkets])
    SumStat[k,3]=std([temp[m][k] for m in 1:nmarkets])
    SumStat[k,4]=minimum([temp[m][k] for m in 1:nmarkets])
    SumStat[k,5]=maximum([temp[m][k] for m in 1:nmarkets])
end
CSV.write(dirtabfig*"/Table11oa.csv",DataFrame(SumStat,:auto))
########################### Analysis of Decomposition
alldata=Matrix(mydata(baseoption,nmarkets));
## location data
uszip=Matrix(combine(groupby(mydata(baseoption,nmarkets),:ziphhn),[:longitude,:latitude] .=> [mean]))

NN=zeros(nmarkets);
for marketid in 1:nmarkets
    alldatatemp=alldata[alldata[:,marketcol].==marketid,[1,2,4]]
    hh=unique(alldatatemp[:,1]); n=length(hh);
    NN[marketid]=n
end
W=NN./sum(NN);

decomdata=Matrix(CSV.read(dirresultsprelim*"/decomp_all_$(nmarkets)_$(default)_1.csv", DataFrame));
nmarkets=Int(maximum(decomdata[:,1])); nsets=Int(maximum(decomdata[:,3]));
Setdist=zeros(nsets,nmarkets); SetdistW=zeros(nsets,nmarkets);

for i in 1:nsets, j in 1:nmarkets
    temp= (decomdata[:,3].==i).*(decomdata[:,1].==j)
    SetdistW[i,j]=decomdata[temp,6][1]*W[j]
    Setdist[i,j]=decomdata[temp,6][1]
end

A=subsets(dY,default); A=vcat(A,sum(SetdistW,dims=2)'); 
A=vcat(A,mean(Setdist,dims=2)'); A=vcat(A,mean(Setdist.>0.0,dims=2)')

######## Number of markets with 5 sets that are considered often
AveM15=mean(sum(Setdist.>0.15,dims=1).<5)
AveM5=mean(sum(Setdist.>0.05,dims=1).<5)

####### Plotting m
Brands=[" CTL ", " GM ", " K ", " O "," Q "]
# Different size
x=1:5
y=[]
for si in 1:5
    y=push!(y,sum(A[6,(sum(A[1:5,:],dims=1)'.==si)[:]]))
end
Figure2oa=bar(x,y, label="", xlabel="|D|", ylabel="Fraction")
savefig(Figure2oa,dirtabfig*"/Figure2oa.pdf")
# |D|=1 and |D|=4
x1=Brands
Sets4=findall((sum(A[1:5,:],dims=1)'.==4)[:])
x4=[prod(Brands[Bool.(A[1:5,Sets4[i]])]) for i in 1:length(Sets4)]
y1=A[6,(sum(A[1:5,:],dims=1)'.==1)[:]]./sum(A[6,(sum(A[1:5,:],dims=1)'.==1)[:]])
y4=A[6,(sum(A[1:5,:],dims=1)'.==4)[:]]./sum(A[6,(sum(A[1:5,:],dims=1)'.==4)[:]])
Figure3oaa=bar(x1,y1, title="|D|=1", label="", xlabel="Sets", ylabel="Fraction")
Figure3oab=bar(x4,y4, title="|D|=4", label="", xlabel="Sets", ylabel="Fraction")
Figure3oa=plot(Figure3oaa,Figure3oab,layout=(2,1))
savefig(Figure3oa,dirtabfig*"/Figure3oa.pdf")

# Brands
x1=Brands
y1=[sum(A[6,findall(A[i,:].==1)]) for i in 1:5]
Figure4oa=bar(x1,y1, label="", xlabel="Brands", ylabel="Fraction")
savefig(Figure4oa,dirtabfig*"/Figure4oa.pdf")
#

# 2 most frequent sets (Online Appendix)
tt=0.05
Actsets=findall(A[6,:].>tt)
x=[prod(Brands[A[1:5,Actsets[i]].==1]) for i in 1:length(Actsets)]
y = A[6,Actsets]


## Plotting F
# Naiv and Estimated shares for full set
#FC
Bigsets=findall(W.>sort(W)[end-5])[1:4]
decomdataNaive=Matrix(CSV.read(dirresultsprelim*"/decomp_all_$(nmarkets)_3_1.csv", DataFrame)); #default=3, which means that the direct shares are computed
FrNaive=zeros(dY,length(Bigsets))
FrEst=zeros(dY,length(Bigsets))

for y in 1:dY, j in 1:length(Bigsets)
    temp= (decomdata[:,2].==1).*(decomdata[:,4].==y).*(decomdata[:,1].==Bigsets[j])
    FrNaive[y,j]=decomdataNaive[temp,5][end]
    FrEst[y,j]=decomdata[temp,5][end]
end

FrNaive=round.(FrNaive,digits=3)
FrEst=round.(FrEst,digits=3)
CSV.write(dirtabfig*"/Table1left.csv",DataFrame(FrNaive,:auto))
CSV.write(dirtabfig*"/Table1right.csv",DataFrame(FrEst,:auto))
# Market 1
FrEstM1=zeros(dY,2)
for y in 1:dY
    temp= (decomdata[:,2].==1).*(decomdata[:,4].==y).*(decomdata[:,1].==1)
    FrEstM1[y,:]=decomdata[temp,5][end-1:end]
end
FrEstM1=round.(FrEstM1,digits=3)
CSV.write(dirtabfig*"/Table12oa.csv",DataFrame(FrEstM1,:auto))

