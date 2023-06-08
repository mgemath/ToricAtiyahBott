export
    a_point


"""
    a_point(v)

The cohomology of a point of `v`.
# Arguments
- `v::NormalToricVariety`: the toric variety.
"""
function a_point(v::NormalToricVariety)::CohomologyClass

    return cohomology_class(rational_equivalence_class(v, [Int64(i == 1) for i in 1:n_cones(v)]))
end

function multinomial(tup::Tuple{Vararg{Int64}})::Int64
    
    SUM::Int64 = 0
    ans::Int64 = 1
    
    @inbounds for i in tup
        SUM += i
        ans *= binomial(SUM, i)
    end

    return ans
end

function stirling_tuple(n::Int64)::Tuple{Vararg{Int64}}
    dp = Array{Int64}(undef, n, n)

    for i = 1:n
        for j = 1:i
            if j == i
                dp[i, j] = 1
            elseif j == 1
                dp[i, j] = factorial(i-1)
            else
                dp[i, j] = dp[i-1, j-1] + (i-1) * dp[i-1, j]
            end
        end
    end

    return (dp[n,:]...,)
end