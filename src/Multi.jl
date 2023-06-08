struct Ms #Multiplicities
    v::NormalToricVariety
    NEF::Vector{CohomologyClass}
    g::Graph{Undirected}
    col::Tuple{Vararg{Int64}}
    inv_curve::Dict{Tuple{Int64,Int64},CohomologyClass}
    beta::CohomologyClass
end

function subtract(M::Ms, marks::Union{Vector{Int64}, Tuple{Vararg{Int64}}})::CohomologyClass
    ans = M.beta

    for (index, e) in enumerate(edges(M.g))
        ans -= (marks[index]*M.inv_curve[M.col[dst(e)],  M.col[src(e)]])
    end

    return ans
end

function multip(M::Ms)::Base.Iterators.Filter

    bound = Vector{Int64}(undef, ne(M.g))
    
    for i in eachindex(bound)
        for j in Iterators.countfrom(1, 1)
            vec = [j*Int(i == k)+Int(i != k) for k in eachindex(bound)]
            if !is_effective_class(M.NEF, subtract(M, vec))
                bound[i] = j-1
                break
            end
        end
    end
    
    return Iterators.filter(m -> istrivial(subtract(M, m)), Iterators.product([1:bound[i] for i in eachindex(bound)]...))
end

function is_effective_class(D::Vector{CohomologyClass}, C::CohomologyClass)::Bool # To be used when D is the vector of nef divisors

    return all(i-> integrate(D[i]*C) > -1, eachindex(D))
end

function nef_divisors(v::NormalToricVariety)::Vector{CohomologyClass}

    ans = CohomologyClass[]
    if dim(v) == 1
        ans = [cohomology_class(v, gens(cohomology_ring(v))[1])]
        return ans
    end

    position = [0 for _ in 1:n_cones(v)]
    len = length(cones(v,1))

    all_curves = CohomologyClass[]
    for index in 1:length(cones(v, dim(v) - 1))
        position[n_maximal_cones(v) + index] = 1
        curve = cohomology_class( rational_equivalence_class(v, position) )
        position[n_maximal_cones(v) + index] = 0
        push!(all_curves, curve)
    end

    for D in Iterators.map(p -> cohomology_class(toric_divisor(v,p)), Iterators.map(i -> [Int(i == j) for j in 1:len], 1:len))
        D in ans && continue

        if all(curve -> integrate(curve*D) > -1, all_curves)
            push!(ans, D)
        end
    end

    return ans
end