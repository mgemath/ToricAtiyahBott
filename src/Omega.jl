function omega(v::NormalToricVariety, n_SIGMA1::Int64, n_SIGMA2::Int64, scalars::Tuple{Vararg{T}})::T

    SIGMA1 = maximal_cones(v)[n_SIGMA1]
    SIGMA2 = maximal_cones(v)[n_SIGMA2]

    ud = QQFieldElem[]

    for ray in rays(SIGMA1)
        ray in rays(SIGMA2) && continue
        for pol_ray in rays(polarize(SIGMA1))
            dot(ray, pol_ray) == 0 && continue
            # ud = pol_ray
            ud = lcm(denominator.(pol_ray)) * pol_ray
            break
        end
        break
    end

    ans = F(0)
    # ans = zero(T)

    for (k, vi) in enumerate(rays(v))# keys(rays_enum)
        # ans += Int64(dot(vi, ud))*scalars[k]
        ans += dot(vi, ud) * scalars[k]
    end

    return ans
end

function omega_dict(v::NormalToricVariety, nc::Dict{Int64,Vector{Int64}}, scalars::Tuple{Vararg{T}})::Dict{Tuple{Int64,Int64},T}

    ans = Dict{Tuple{Int64,Int64},T}()

    for sigma1 in keys(nc)
        for sigma2 in nc[sigma1]
            ans[(sigma1, sigma2)] = omega(v, sigma1, sigma2, scalars)
        end
    end

    return ans
end

function omega_total(od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, sigma::Int64)::T

    return prod(gamma -> od[(sigma, gamma)], nc[sigma])
end

function omega_F(od::Dict{Tuple{Int64,Int64},T}, sigma1::Int64, sigma2::Int64, d::Int64)::T

    return od[(sigma1, sigma2)] // d
end