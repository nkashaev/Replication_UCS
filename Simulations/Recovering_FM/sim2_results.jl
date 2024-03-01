#This files generates Tables 3-10 in the online appendix
using Pkg
Pkg.activate(".")
using CSV, DataFrames
using LinearAlgebra
########################### Dir ################################################
tempdir1=@__DIR__
rootdir=tempdir1[1:findfirst("Replication_UCS",tempdir1)[end]]
dir=rootdir*"/Simulations/Recovering_FM/"
dirresults=dir*"Results_sim2/"
########################### Functions #########################################
function results_sim2(Pd,Pyd,model,sampsize)
    numchoices,numsets=size(Pyd)
    M=Matrix(CSV.read(dirresults*"M_$(model)_$(sampsize).csv", DataFrame))
    F=Matrix(CSV.read(dirresults*"F_$(model)_$(sampsize).csv", DataFrame))
    Mbias=zeros(length(M[:,1]))
    Mrmse=zeros(length(M[:,1]))
    Fbias=zeros(length(F[:,1]))
    Frmse=zeros(length(F[:,1]))
    nsim=size(M,2)
    for s in 1:nsim
      Mbias=Mbias .+ (M[:,s]-Pd)./nsim
      Fbias=Fbias .+ (F[:,s]-vcat(Pyd[:],Pyd[:],Pyd[:]))./nsim
      Mrmse=Mrmse .+ (M[:,s]-Pd).^2 ./nsim
      Frmse=Frmse .+ (F[:,s]-vcat(Pyd[:],Pyd[:],Pyd[:])).^2 ./nsim
    end
  
    return Mbias,Fbias,sqrt.(Mrmse),sqrt.(Frmse)
end
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
    
#################################################################
SS=[2000, 5000, 10000, 50000]
dn=length(SS)
BiasMDgp1=zeros(dn,5)
RmseMDgp1=zeros(dn,5)
BiasMDgp2=zeros(dn,5)
RmseMDgp2=zeros(dn,5)

BiasFDgp1=zeros(dn,4)
RmseFDgp1=zeros(dn,4)
BiasFDgp2=zeros(dn,10)
RmseFDgp2=zeros(dn,10)

model="core"
Pd,Pyd=dgpparam(model)
Rep=(Pyd.>0.0)
Rep[1,:]=zeros(5)
for i in 1:dn
    Mbias,Fbias,Mrmse,Frmse=results_sim2(Pd,Pyd,model,SS[i])
    BiasMDgp1[i,:]=Mbias'
    RmseMDgp1[i,:]=Mrmse'
    BiasFDgp1[i,:]=reshape(Fbias[1:25],5,5)[Rep]'
    RmseFDgp1[i,:]=reshape(Frmse[1:25],5,5)[Rep]'
end

model="nested"
Pd,Pyd=dgpparam(model)
Rep=(Pyd.>0.0)
Rep[1,:]=zeros(5)
for i in 1:dn
    Mbias,Fbias,Mrmse,Frmse=results_sim2(Pd,Pyd,model,SS[i])
    BiasMDgp2[i,:]=Mbias'
    RmseMDgp2[i,:]=Mrmse'
    BiasFDgp2[i,:]=reshape(Fbias[1:25],5,5)[Rep]'
    RmseFDgp2[i,:]=reshape(Frmse[1:25],5,5)[Rep]'
end

Table3=DataFrame(10^5 * BiasMDgp1,:auto);
Table4=DataFrame(10^3 * RmseMDgp1,:auto);
Table5=DataFrame(10^3 * BiasMDgp2,:auto);
Table6=DataFrame(10^2 * RmseMDgp2,:auto);
Table7=DataFrame(10^4 * BiasFDgp1,:auto);
Table8=DataFrame(10^2 * RmseFDgp1,:auto);
Table9=DataFrame(10^2 * BiasFDgp2,:auto);
Table10=DataFrame(10^2 * RmseFDgp2,:auto);

CSV.write(dir*"Tables/table3.csv", round.(Table3,digits=1))
CSV.write(dir*"Tables/table4.csv", round.(Table4,digits=1))
CSV.write(dir*"Tables/table5.csv", round.(Table5,digits=1))
CSV.write(dir*"Tables/table6.csv", round.(Table6,digits=1))
CSV.write(dir*"Tables/table7.csv", round.(Table7,digits=1))
CSV.write(dir*"Tables/table8.csv", round.(Table8,digits=1))
CSV.write(dir*"Tables/table9.csv", round.(Table9,digits=1))
CSV.write(dir*"Tables/table10.csv", round.(Table10,digits=1))

