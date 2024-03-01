# Functions used in the 2nd set of simulations
# Generating an output of one simulation
function oneMC(seed)
  # Generating data
  data=dgp1(sampsize, Pyd, Pd, seed)
  # Computing the 3-d array of joint probabilites of 3 outcomes
  Pijk=Pijk_sets(data)
  # Final result
  F1s2, F2s2, F3s2, Ms2=step2c(Pijk, Pyd, Pd, Pyd.>0)

  return F1s2, F2s2, F3s2, Ms2
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


# This function computes the final estimators of Frum and M
# Pijk -- 3-d array of joint probabilites of 3 outcomes
# Sets -- the recovered sets
function step2c(Pijk, Pyd, Pd, Sets)
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
  for i=1:Tc
    set_start_value(m[i],Pd[i])
    for j=1:T
      set_start_value(f1[j,i], Pyd[j,i])
      set_start_value(f2[j,i], Pyd[j,i])
      set_start_value(f3[j,i], Pyd[j,i])
    end
  end
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

