# Functions used in the 1 set of simulations

# Generating the output of one simulation
function oneMC(seed)
  # Generating data
  #println(seed)
  data=dgp1(sampsize, Pyd, Pd, seed)
  # Computing the 3-d array of joint probabilites of 3 outcomes
  Pijk=Pijk_sets(data)
  # Estimated number of choice sets
  #d=rankest(reshape(sum(Pijk,dims=3),size(Pijk,1),size(Pijk,1)),1/size(data,1))
  d=rank(reshape(sum(Pijk,dims=3),size(Pijk,1),size(Pijk,1)))
  # Step 1 estimator
  F1s1, F2s1, F3s1, Ms1=step1(Pijk, d, tau)
  # All possible choice sets
  A=subsets(size(Pijk,1),1)
  # Unconstrained by the number of choice sets estimation
  F1, F2, F3 = step2a(Pijk, A, F1s1, F2s1, F3s1, Ms1)
  # Mixed-integer optimization to recover sets
  b=step2b(Pijk,F1,F2,F3,d)
  # Final result
  F1s2s, F2s2s, F3s2s, Ms2s=step2c(Pijk, A[:,(round.(Int,b).==1)])
  # Transforming into large matrices
  F1s2=zeros(size(A)); F2s2=zeros(size(A)); F3s2=zeros(size(A));
  F1s2[1,:]=ones(size(A,2)); F2s2[1,:]=ones(size(A,2)); F3s2[1,:]=ones(size(A,2));
  Ms2=zeros(size(A,2))
  for i in 1:size(A,2)
    for i1 in 1:size(F1s2s,2)
      if A[:,i] == (F1s2s[:,i1].>0.0)
        F1s2[:,i]=F1s2s[:,i1]
        F2s2[:,i]=F2s2s[:,i1]
        F3s2[:,i]=F3s2s[:,i1]
        Ms2[i]=Ms2s[i1]
      end
    end
  end

  return F1s1, F2s1, F3s1, Ms1, F1s2, F2s2, F3s2, Ms2
end

# This function generates data
# sampsize -- number of observations
# Pyd -- matrix of F^{RUM}: Pyd[y,D]=Pr(y|D)
function dgp1(sampsize, Pyd, Pd, seed)
  rng1 = MersenneTwister(seed); rng2 = MersenneTwister(seed^2);
  datagen=zeros(sampsize,3)
  D=rand(rng1,Categorical(Pd),sampsize)
  for i in 1:sampsize
    datagen[i,:]=rand(rng2,Categorical(Pyd[:,D[i]]),3)'
  end
  return datagen
end

# This function computs the 3-d array of joint probabilites of 3 outcomes
# from data
function Pijk_sets(data)
    N=size(data,1)
    dY=Int(maximum(data))
    Pijk=zeros(dY,dY,dY)
    for i in 1:dY, j in 1:dY, k in 1:dY
        Pijk[i,j,k]= mean((data[:,1].==i).*(data[:,2].==j).*(data[:,3].==k))
    end
    return Pijk
end


# This function computes Step-1 estimators of m and Frum.
# Pijk -- 3-d array of joint probabilites of 3 outcomes
# d -- estimated number of choice sets
# taun -- trimming parameter
function step1(Pijk, d, taun)
  f1, f2, f3, m=step1_optimization(Pijk,d)
  T=size(f1,1)
  SSets=[((f1.>taun).*(f2.>taun).*(f3.>taun))[:,i] for i in 1:d]
  Aall=subsets(T,1)
  F1=zeros(size(Aall)); F2=zeros(size(Aall)); F3=zeros(size(Aall))
  F1[1,:]=ones(size(Aall,2)); F2[1,:]=ones(size(Aall,2)); F3[1,:]=ones(size(Aall,2))
  M=zeros(size(Aall,2))
  for i in 1:size(Aall,2)
    for i1 in 1:length(SSets)
      Stemp=vcat(1, SSets[i1][2:end])
      if Aall[:,i] == Stemp
        F1[:,i]=f1[:,i1].*Stemp
        F2[:,i]=f2[:,i1].*Stemp
        F3[:,i]=f3[:,i1].*Stemp
        M[i]=m[i1].*(m[i1].>taun)
      end
    end
  end
  return  F1./sum(F1,dims=1), F2./sum(F2,dims=1), F3./sum(F3,dims=1), M./sum(M)
end

# This function does optimization for Step-1 estimators of m and Frum.
# Pijk -- 3-d array of joint probabilites of 3 outcomes
# d -- estimated number of choice sets
function step1_optimization(Pijk,d)
  T=size(Pijk,1)
  Tc=d
  mstep1 = Model(KNITRO.Optimizer)
  set_optimizer_attribute(mstep1,"outlev",0)
  set_optimizer_attribute(mstep1,"ms_enable",1)
  @variable(mstep1, 0 <= m[1:Tc] <=1)
  @variable(mstep1, 0 <= f1[1:T,1:Tc] <= 1)
  @variable(mstep1, 0 <= f2[1:T,1:Tc] <= 1)
  @variable(mstep1, 0 <= f3[1:T,1:Tc] <= 1)
  @constraint(mstep1, addm,sum(m[t] for t in 1:Tc)==1)
  @constraint(mstep1, addf1[t in 1:Tc],sum(f1[j,t] for j in 1:T)==1)
  @constraint(mstep1, addf2[t in 1:Tc],sum(f2[j,t] for j in 1:T)==1)
  @constraint(mstep1, addf3[t in 1:Tc],sum(f3[j,t] for j in 1:T)==1)
  @NLobjective(mstep1, Max, -sum(sum(sum((Pijk[i,j,k]-sum(f1[i,t]*f2[j,t]*f3[k,t]*m[t] for t in 1:Tc))^2 for i in 1:T) for j in 1:T) for k in 1:T))
  JuMP.optimize!(mstep1)
  return value.(f1), value.(f2), value.(f3), value.(m)
end

# This function does first optimization for Step-2 estimators of m and Frum.
# Pijk -- 3-d array of joint probabilites of 3 outcomes
# A is the matrix of all possible choice sets
# Fts1, Ms1 -- Step-1 estimators
function step2a(Pijk, A, F1s1, F2s1, F3s1, Ms1)
  T,Tc=size(A)
  mstep2a = Model(KNITRO.Optimizer)
  set_optimizer_attribute(mstep2a,"outlev",0)
  @variable(mstep2a, 0 <= m[1:Tc] <=1)
  @variable(mstep2a, 0 <= f1[1:T,1:Tc] <= 1)
  @variable(mstep2a, 0 <= f2[1:T,1:Tc] <= 1)
  @variable(mstep2a, 0 <= f3[1:T,1:Tc] <= 1)
  @constraint(mstep2a, addm,sum(m[t] for t in 1:Tc)==1)
  @constraint(mstep2a, addf1[t in 1:Tc],sum(f1[j,t] for j in 1:T)==1)
  @constraint(mstep2a, addf2[t in 1:Tc],sum(f2[j,t] for j in 1:T)==1)
  @constraint(mstep2a, addf3[t in 1:Tc],sum(f3[j,t] for j in 1:T)==1)
  @constraint(mstep2a, zerof1[i in 1:T, t in 1:Tc], f1[i,t] <= A[i,t])
  @constraint(mstep2a, zerof2[i in 1:T, t in 1:Tc], f2[i,t] <= A[i,t])
  @constraint(mstep2a, zerof3[i in 1:T, t in 1:Tc], f3[i,t] <= A[i,t])
  for i=1:Tc
    set_start_value(m[i],Ms1[i])
    for j=1:T
      set_start_value(f1[j,i], F1s1[j,i])
      set_start_value(f2[j,i], F2s1[j,i])
      set_start_value(f3[j,i], F3s1[j,i])
    end
  end
  @NLobjective(mstep2a, Max, -sum(sum(sum((Pijk[i,j,k]-sum(f1[i,t]*f2[j,t]*f3[k,t]*m[t] for t in 1:Tc))^2 for i in 1:T) for j in 1:T) for k in 1:T))
  JuMP.optimize!(mstep2a)
  return value.(f1), value.(f2), value.(f3) #objective_value(mstep2a)
end

# This function recovers sets
# Pijk -- 3-d array of joint probabilites of 3 outcomes
# Ft -- preliminary Step 2 estimators of Frum
# d -- estimated number of choice sets
function step2b(Pijk,F1,F2,F3,d)
  T,Tc=size(F1)
  B=zeros(T,T,T,Tc)
  for i in 1:T, j in 1:T, k in 1:T, cs in 1:Tc
    B[i,j,k,cs]=F1[i,cs]*F2[j,cs]*F3[k,cs]
  end
  mstep2b = Model(KNITRO.Optimizer)
  set_optimizer_attribute(mstep2b,"outlev",0)
  @variable(mstep2b, 0 <= m[1:Tc] <=1)
  @variable(mstep2b, y[i=1:Tc], Bin)
  @constraint(mstep2b, addm,sum(m[t] for t in 1:Tc)==1)
  @constraint(mstep2b, addy,sum(y[t] for t in 1:Tc)<=d)
  @constraint(mstep2b, zerom[t in 1:Tc], m[t] <=y[t])
  @objective(mstep2b, Max, -sum(sum(sum((Pijk[i,j,k]-sum(B[i,j,k,t]*m[t] for t in 1:Tc))^2 for i in 1:T) for j in 1:T) for k in 1:T))
  JuMP.optimize!(mstep2b)
  return value.(y)
end

# This function computes the final estimators of Frum and M
# Pijk -- 3-d array of joint probabilites of 3 outcomes
# Sets -- the recovered sets
function step2c(Pijk, Sets)
  T,Tc=size(Sets)
  mstep2c = Model(KNITRO.Optimizer)
  set_optimizer_attribute(mstep2c,"outlev",0)
  @variable(mstep2c, 0 <= m[1:Tc] <=1)
  @variable(mstep2c, 0 <= f1[1:T,1:Tc] <= 1)
  @variable(mstep2c, 0 <= f2[1:T,1:Tc] <= 1)
  @variable(mstep2c, 0 <= f3[1:T,1:Tc] <= 1)
  @constraint(mstep2c, addm, sum(m[t] for t in 1:Tc)==1)
  @constraint(mstep2c, addf1[t in 1:Tc], sum(f1[j,t] for j in 1:T)==1)
  @constraint(mstep2c, addf2[t in 1:Tc], sum(f2[j,t] for j in 1:T)==1)
  @constraint(mstep2c, addf3[t in 1:Tc], sum(f3[j,t] for j in 1:T)==1)
  @constraint(mstep2c, zerof1[i in 1:T, t in 1:Tc], f1[i,t] <=Sets[i,t])
  @constraint(mstep2c, zerof2[i in 1:T, t in 1:Tc], f2[i,t] <=Sets[i,t])
  @constraint(mstep2c, zerof3[i in 1:T, t in 1:Tc], f3[i,t] <=Sets[i,t])
  @NLobjective(mstep2c, Max, -sum(sum(sum((Pijk[i,j,k]-sum(f1[i,t]*f2[j,t]*f3[k,t]*m[t] for t in 1:Tc))^2 for i in 1:T) for j in 1:T) for k in 1:T))
  JuMP.optimize!(mstep2c)
  return value.(f1), value.(f2), value.(f3), value.(m)
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

