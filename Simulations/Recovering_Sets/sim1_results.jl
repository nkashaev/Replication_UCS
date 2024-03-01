#This files generates Tables 1-2 in the online appendix
using Pkg
Pkg.activate(".")
using CSV, DataFrames
using LinearAlgebra
########################### Dir ################################################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Replication_UCS",tempdir1)[end]]
dir=rootdir*"/Simulations/Recovering_Sets/"
dirresults=dir*"Results_sim1/"
########################### Functions #########################################
function dgpparam(model)
    ## Distribution over the sets
    Pd=[0.2, 0.15, 0.3, 0.15, 0.2]
    #Distribution over options conditional on sets:
    #Pyd_number of options_number of sets
    Pyd_nested_5_5=[1.00 0.60 0.50 0.25 0.10;
                    0.00 0.40 0.20 0.35 0.25;
                    0.00 0.00 0.30 0.25 0.15;
                    0.00 0.00 0.00 0.15 0.30;
                    0.00 0.00 0.00 0.00 0.20]
  
    Pyd_core_5_5=[1.00 0.60 0.50 0.40 0.20;
                  0.00 0.40 0.00 0.00 0.00;
                  0.00 0.00 0.50 0.00 0.00;
                  0.00 0.00 0.00 0.60 0.00;
                  0.00 0.00 0.00 0.00 0.80]
    ## Model
      if model=="nested"
        Pyd=Pyd_nested_5_5
      elseif model=="core"
        Pyd=Pyd_core_5_5
      end
    return Pd, Pyd
end
# This function computes all nonempty subsets of a choice set of size dP that all contain default
# If Default=0, then there is no defualt
function subsets(dP,default)
    default==0 ? dA=dP : dA=dP-1
    A=zeros(Int8,2^dA,dA)
    indexA=ones(Int8,dA,1)
    p=1
    ready = false
    while ~ready
      A[p,:]=indexA'
      p=p+1
      ready = true
      for a = dA:-1:1
        indexA[a] = indexA[a] + 1
        if indexA[a] <= 2
          ready = false
          break
        end
        indexA[a] = 1
      end
    end
    A=(A.-1)'
    default==0 ? A=A[:,2:end] : A=vcat(A[1:default-1,:],ones(Int,size(A,2))',A[default:end,:])
    return A
end
  
numsets=5
numchoices=5  
function results_sim1(Pd,Pyd,model,sampsize)
    numchoices,numsets=size(Pyd)
    Bs1=Matrix(CSV.read(dirresults*"Ms1_$(model)_$(sampsize).csv", DataFrame)).> 0.0
    Bs2=Matrix(CSV.read(dirresults*"Ms2_$(model)_$(sampsize).csv", DataFrame)).> 0.0
    Aall=subsets(numchoices,1)
    btrue=zeros(size(Aall,2))
    for i in 1:length(btrue), j in 1:length(Pd)
      if Aall[:,i]==(Pyd[:,j].>0.0)
        btrue[i]=1.0
      end
    end
    atrue=findall(btrue.==1)
    inter1=0.0;inter2=0.0
    s1=0.0; s2=0.0;
    s1dist =zeros(numsets); s2dist = zeros(numsets);
    nsim=size(Bs1,2)
    for s in 1:nsim
      s1=s1+(btrue==Bs1[:,s])/nsim
      s2=s2+(btrue==Bs2[:,s])/nsim
      inter1=inter1+length(intersect(atrue,findall(Bs1[:,s].==1)))/nsim
      inter2=inter2+length(intersect(atrue,findall(Bs2[:,s].==1)))/nsim
    end
    s1dist=1.0 .- sum(Bs1[atrue,:],dims=2)./nsim
    s2dist=1.0 .- sum(Bs2[atrue,:],dims=2)./nsim
    return round(s1,digits=3),round(s2,digits=3),round(inter1,digits=3),round(inter2,digits=3),round.(s1dist,digits=3),round.(s2dist,digits=3)
end
  
##################################################
SS=[2000,5000,10000,50000]
numsets=5
numchoices=5

Table1=zeros(4,4)
Table2=zeros(4,4)

model="core"
Pd,Pyd=dgpparam(model)
for i in 1:length(SS)
    s1,s2,iner1,iner2,s1dist,s2dist=results_sim1(Pd,Pyd,model,SS[i])
    Table1[1,i]=s1*100
    Table1[2,i]=s2*100
    Table2[1,i]=iner1
    Table2[2,i]=iner2
end

model="nested"
Pd,Pyd=dgpparam(model)
for i in 1:length(SS)
    s1,s2,iner1,iner2,s1dist,s2dist=results_sim1(Pd,Pyd,model,SS[i])
    Table1[3,i]=s1*100
    Table1[4,i]=s2*100
    Table2[3,i]=iner1
    Table2[4,i]=iner2
end
CSV.write(dir*"Tables/table1.csv", DataFrame(Table1,:auto))
CSV.write(dir*"Tables/table2.csv", round.(DataFrame(Table2,:auto),digits=2))