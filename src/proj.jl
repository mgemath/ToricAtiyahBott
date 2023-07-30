export
    proj
    
    
"""
    proj(E...)

Projectivization of the direct sum of the line bundles (or toric divisors) in E. 
# Arguments
- `E`: a sequence of line bundles (or toric divisors) on a toric variety.


!!! note

    The line bundles must be on the same toric variety.


# Example
Let us construct ``v = \\mathbb{P}(\\mathcal{O}_{\\mathbb{P}^1\\times\\mathbb{P}^1}(1,0)\\oplus\\mathcal{O}_{\\mathbb{P}^1\\times\\mathbb{P}^1}(0,1))``.
```julia-repl
julia> P1 = projective_space(NormalToricVariety, 1);
julia> D0 = toric_divisor(P1, [0,0]);
julia> P1xP1 = proj(D0, D0);
julia> O10 = toric_line_bundle(P1xP1, [1,0]);
julia> O01 = toric_line_bundle(P1xP1, [0,1]);
julia> v = proj(O10, O01);
```
"""
function proj(E...)::Union{Nothing, NormalToricVariety}

    v = toric_variety(E[1])

    length(E) == 1 && return v

    if any(i -> toric_variety(E[i]) != v, eachindex(E))
        error("The divisors are defined on different toric varieties.")
        return nothing
    end

    PF_Pr = normal_fan(simplex(length(E) - 1))
    l = rays(PF_Pr)

    modified_ray_gens = Dict{RayVector{QQFieldElem}, RayVector{QQFieldElem}}()

    for sigma in maximal_cones(v)
        for ray in rays(sigma)
            ray in keys(modified_ray_gens) && continue
            modified_ray_gens[ray] = vcat(ray, -sum(i -> dot(m_sigma(sigma, E[i]),ray)*l[i], eachindex(E)))
        end
        length(keys(modified_ray_gens)) == nrays(v) && break
    end

    # new_maximal_cones = Vector{Int64}[]
    new_maximal_cones = Vector{Vector{Int64}}(undef, n_maximal_cones(v)*length(E))
    index = 1

    for a in 1:n_maximal_cones(v)
        first = [row(ray_indices(maximal_cones(v)), a)...,]
        for b in eachindex(E)
            second = [row(ray_indices(maximal_cones(PF_Pr)), b)...,] .+ nrays(v)
            # push!(new_maximal_cones, vcat(first, second))
            new_maximal_cones[index] = vcat(first, second)
            index += 1
        end
    end

    total_rays_gens = vcat([modified_ray_gens[ray] for ray in rays(v)], [vcat(RayVector(zeros(Int64, dim(v))), l[i]) for i in eachindex(E)])

    return NormalToricVariety(polyhedral_fan(total_rays_gens, IncidenceMatrix(new_maximal_cones)))
end

function m_sigma(sigma::Cone{QQFieldElem}, D::Union{ToricDivisor,ToricLineBundle})::RayVector{QQFieldElem}
    
    ans = RayVector(zeros(QQFieldElem, dim(sigma)))
    coeff = coefficients(isa(D,ToricDivisor) ? D : toric_divisor(D))

    rays_pol = rays(polarize(sigma))
    dual_ray = QQFieldElem[]

    for ray in rays(sigma)
        for pol_ray in rays_pol
            dot(ray, pol_ray) == 0 && continue
            dual_ray = lcm(denominator.(pol_ray))*pol_ray
            break
        end
        i = index_in(ray, rays(toric_variety(D)))
        ans += -coeff[i]*dual_ray
    end

    return ans
end