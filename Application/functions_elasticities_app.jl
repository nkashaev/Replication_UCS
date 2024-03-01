# This function computes GMM estimator of the price coefficient with identity weight matrix
# data          -- DataFrame
# baseprod      -- the base product
# mstart        -- set to 1 to have multistart, to 0 otherwise
function gmmbeta(data, baseprod, mstart=1)
    maxnopt=maximum(data.product_ids)                       # Total number of alternatives
    choiceset=findall(data.shares[1:maxnopt].>0.0)          # Choice set
    choicesetnobase=setdiff(choiceset,baseprod)             # Choice set minus the base product
    Nchoicesetnobase=length(choicesetnobase)                # Cardinality of the choice set without the base product
    if Nchoicesetnobase==length(choiceset)
        return "Base option does not belong to the choice set"
    end
    def=data.product_ids .== baseprod                       # Observations that correspond to the base product

    # Setting up the optimization
    gmmopt = Model(KNITRO.Optimizer)
    set_optimizer_attribute(gmmopt,"outlev",0)              # Turning off the ouput
    set_optimizer_attribute(gmmopt,"ms_enable",mstart)      # Multistart option
    @variable(gmmopt, alpha[1:4,1:Nchoicesetnobase])        # Coefficients in front of 4 controls: constant, age, income, hhsize
    @variable(gmmopt, beta)                                 # Price coefficient
    # Setting up dependent variable (y), relative price (p),
    # instruments (Z),  and controls (Xnp)
    y=zeros(sum(def),maxnopt)
    p=zeros(sum(def),maxnopt)
    Z=zeros(sum(def),size(alpha,1)+2,maxnopt)               # Instruments include controls + 2 excluded instruments
    Xnp=zeros(sum(def),size(alpha,1),maxnopt)
    for prod in choicesetnobase
        pr=(data.product_ids.==prod)                        # Observations that correspond to prod
        y[:,prod]=log.(data.shares[pr]./data.shares[def])   # Log shares
        p[:,prod]=data.price[pr].-data.price[def]           # Relative prices
        # Controls inculde constant, age, income, hhsize
        Xnp[:,:,prod]=[ones(sum(pr)) Matrix(data[pr,11:13])]
        # Instruments contain controls,
        # relative average sizes of compiting products within a market,
        # and average prices of the same product in different markets
        Z[:,:,prod]=[ones(sum(pr)) data.demand_instruments1[pr] .- data.demand_instruments1[def] data.demand_instruments0[pr] .- data.demand_instruments0[def] Matrix(data[pr,11:13])]
    end
    # Setting outliers to be zero
    for prod in choicesetnobase
        outl=(abs.(y[:,prod]).>2.0)                           # The threshold is 2.0
        y[outl,prod].=0.0
        p[outl,prod].=0.0
        Z[outl,:,prod].=0.0
        Xnp[outl,:,prod].=0.0
    end
    ninst=size(Z,2)                                         # Number of instruments
    # Moments
    moment=@expression(gmmopt, [sum( Z[i,m,prod]*(y[i,prod] - beta*p[i,prod] - sum(Xnp[i,k,prod]*alpha[k,findall(choicesetnobase.==prod)[1]] for k=1:size(alpha,1))) for i=1:sum(def))/sum(def) for m=1:ninst, prod in choicesetnobase][:])
    # Objective
    @objective(gmmopt, Min, moment'*moment)
    JuMP.optimize!(gmmopt)
    b1=value.(beta); a1=value.(alpha);
    moment1=[sum( Z[i,m,prod]*(y[i,prod] - b1*p[i,prod] - sum(Xnp[i,k,prod]*a1[k,findall(choicesetnobase.==prod)[1]] for k=1:size(alpha,1))) for i=1:sum(def))/sum(def) for m=1:ninst, prod in choicesetnobase][:]
    #Variance of the moments and Jacobian
    Vm=zeros(ninst*Nchoicesetnobase,ninst*Nchoicesetnobase)
    G=zeros(ninst*Nchoicesetnobase,prod(size(a1))+1)
    for i in 1:sum(def)
        mi=Z[i,:,choicesetnobase][:].*repeat(y[i,choicesetnobase] .- b1*p[i,choicesetnobase] .- [Xnp[i,:,prod]'*a1[:,findall(choicesetnobase.==prod)[1]] for prod in choicesetnobase],1,ninst)'[:]
        Vm=Vm + mi*mi'./sum(def)
        cc1=zeros(ninst*Nchoicesetnobase,prod(size(a1)))
        cc2=zeros(ninst*Nchoicesetnobase)
        k=0
        for  pp in choicesetnobase
            zi=Z[i,:,pp]
            xi=Xnp[i,:,pp]
            cc1[k*length(zi)+1:(k+1)*length(zi),k*length(xi)+1:(k+1)*length(xi)]=zi*xi'
            cc2[k*length(zi)+1:(k+1)*length(zi),:]=zi*p[i,pp]
            k=k+1
        end
        G=G+hcat(cc2,cc1)./sum(def)
    end
    Vm=Vm-moment1*moment1'
    Weight=inv(Vm)
    @objective(gmmopt, Min, moment'*Weight*moment)
    JuMP.optimize!(gmmopt)

    b1=value.(beta); a1=value.(alpha);
    moment1=[sum( Z[i,m,prod]*(y[i,prod] - b1*p[i,prod] - sum(Xnp[i,k,prod]*a1[k,findall(choicesetnobase.==prod)[1]] for k=1:size(alpha,1))) for i=1:sum(def))/sum(def) for m=1:ninst, prod in choicesetnobase][:]
    #Variance of the moments and Jacobian
    Vm=zeros(ninst*Nchoicesetnobase,ninst*Nchoicesetnobase)
    for i in 1:sum(def)
        mi=Z[i,:,choicesetnobase][:].*repeat(y[i,choicesetnobase] .- b1*p[i,choicesetnobase] .- [Xnp[i,:,prod]'*a1[:,findall(choicesetnobase.==prod)[1]] for prod in choicesetnobase],1,ninst)'[:]
        Vm=Vm + mi*mi'./sum(def)
    end
    Vm=Vm-moment1*moment1'
    se=sqrt(Complex((inv(G'*inv(Vm)*G)[1,1]/sum(def))))
    if real(se)==0.0
        display("error in se due to noninvertibility")
    end
    return b1, real(se)
end


#This function generates the dataframe needed for gmmbeta function
function data_gen(set_number, baseoption, nmarkets, default, distrank)
    alldata=Matrix(mydata(baseoption, nmarkets))
    T=Int(maximum(alldata[:,2])); dY=Int(maximum(alldata[:,4]));
    dY==5 ? marketcol=15 : marketcol=13
    nmarkets=length(unique(alldata[:,marketcol]))
    #Prices rescaled by the size
    pricesall=zeros(T,nmarkets,dY)
    for t in 1:T, marketid in 1:nmarkets, prod in 1:dY
        ind=(alldata[:,2].==t).*(alldata[:,4].==prod).*(alldata[:,marketcol].==marketid)
        pricesall[t,marketid,:]=unique(alldata[ind,marketcol-2dY:marketcol-2dY+dY-1],dims=1)./unique(alldata[ind,marketcol-2dY+dY:marketcol-2dY+dY+dY-1],dims=1)
    end

    #Finding geografical centers of markets
    coordinates=zeros(nmarkets,2)
    for i in 1: nmarkets
        coordinates[i,:]=[mean(unique(alldata[alldata[:,marketcol].==i,marketcol+1])) mean(unique(alldata[alldata[:,marketcol].==i,marketcol+2]))]
    end

    #Neighbours
    neighbours=[[] for i in 1:nmarkets]
    for i in 1:nmarkets
        Distij=[norm(coordinates[i,:] .- coordinates[j,:]) for j in 1:nmarkets]
        neighbours[i]=setdiff(findall(Distij.<=sort(Distij)[distrank+1]),[i])
    end

    #Decomposed data
    decomdata=Matrix(CSV.read(dirresultsprelim*"/decomp_all_$(nmarkets)_$(default)_$(baseoption).csv", DataFrame))
    decomdata=decomdata[decomdata[:,6].>0.0 ,:] # Dropping inactive sets
    decomdata=decomdata[decomdata[:,3].==set_number,:] # Picking only set_number
    decomdata=unique(decomdata,dims=1)
    Markets=Int.(unique(decomdata[:,1]))

    Output=DataFrame(market_ids=Float64[], trip_ids=Int[], set_ids=Int[], product_ids=Int[], shares=Float64[], set_pr=Float64[], price=Float64[], demand_instruments0=Float64[], demand_instruments1=Float64[], size=Float64[], age=Float64[], income=Float64[], hhsize=Float64[])
    for i in 1:size(decomdata,1)
        t=Int(decomdata[i,2])
        market=Int(decomdata[i,1])
        prod=Int(decomdata[i,4])
        Ind=(alldata[:,2].==t).*(alldata[:,marketcol].==market)
        c=[decomdata[i,:]' pricesall[t,market,prod] mean(pricesall[t,neighbours[market],prod]) mean(alldata[Ind,marketcol-6 .+ setdiff(collect(1:5),prod)]) mean(alldata[Ind,marketcol-6+prod]) mean(alldata[Ind,marketcol+3]) mean(alldata[Ind,marketcol+4]) mean(alldata[Ind,marketcol+5])]
        push!(Output,c)
    end

    return Output
end



# This function computes all nonempty subsets of a choice set of size dP that all contain alternative 1
# If default=0, then there is no option that is always considered
# If default=1, then alternative 1 is always considered
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
    default==0 ? A=A[:,2:end] : A=vcat(ones(Int,size(A,2))',A)
    return A
end

## This function reads the data and puts the baseoption as the first one and keeps the order of the rest the same
# baseoption -- this option will be put first
# nmarkets   -- number of markets determined by K-means
function mydata(baseoption, nmarkets)
    alldataframe=CSV.read(dirdata*"/extract_2016_2018_$(nmarkets).csv", DataFrame) # Original data
    #Turning the baseoption to 1, k<baseoption to k+1, keepeing k if k>baseoption
    for i in 1:length(alldataframe.prod)
        if alldataframe.prod[i]==baseoption
            alldataframe.prod[i]=1
        elseif alldataframe.prod[i]<baseoption
            alldataframe.prod[i]=alldataframe.prod[i]+1
        end
    end
    #Shifting the baseoption price and size to the first positions
    return alldataframe[!, vcat(collect(1:4), 4+baseoption,4 .+ setdiff(collect(1:5), baseoption),
                         9+baseoption, 9 .+ setdiff(collect(1:5), baseoption), collect(15:20))]
end
