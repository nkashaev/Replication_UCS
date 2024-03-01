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
