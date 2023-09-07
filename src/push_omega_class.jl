export
    push_omega


"""
    push_omega(l)

Equivariant class of the push-forward under the forgetful map of ``ev^*l`` tensored the cotangent bundle of the forgetful map.
# Arguments
- `l::ToricLineBundle`: the line bundle.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{3},1)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{2}\\cdot\\mathrm{ev}_{2}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{3}\\cdot\\mathrm{c_{top}}(\\delta_{*}(\\omega_{\\delta}\\otimes\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(2))) &= 1 \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using Oscar, ToricAtiyahBott)
julia> v = projective_space(NormalToricVariety, 3);

julia> line = cohomology_class(toric_line_bundle(v, [1]))^2;

julia> P = push_omega(toric_line_bundle(v, [2]))*ev(1, line)*ev(2, a_point(v));

julia> IntegrateAB(v, line, 2, P, show_bar=false);
Result: 1
```
"""
function push_omega(l::ToricLineBundle)::EquivariantClass

    cc = cohomology_class(toric_divisor(l))
    # SPLIT = split_cc(cc)
    # Z = [(i, Int64.(exponents(SPLIT[i])[1,:]), Int64(Oscar.coefficients(SPLIT[i])[1])) for i in eachindex(SPLIT)]

    rule = :(_push_omega(v, od, nc, iv, g, col, weights, marks, $cc))
    # rule = :(_push_ev(v, od, nc, iv, g, col, weights, marks, $l))

    return EquivariantClass(rule, eval(:((v, od, nc, iv, g, col, weights, marks) -> $rule)))
end

# function _push_ev(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, cc::CohomologyClass, Z::Vector{Tuple{Int64, Matrix{Int64}, Int64}})::T
function _push_omega(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, cc::CohomologyClass)::T
    ans = F(1)

    # evd = Dict{Int64, T}()

    for vert in 1:nv(g)
        if !((col[vert], cc) in keys(v.__attrs[:ev_dict]))
            v.__attrs[:ev_dict][(col[vert], cc)] = just_ev(v, od, nc, col[vert], cc)
        end

        1 - length(all_neighbors(g, vert)) == 0 && continue
        v.__attrs[:ev_dict][(col[vert], cc)] == 0 && return v.__attrs[:ev_dict][(col[vert], cc)]
        ans *= (v.__attrs[:ev_dict][(col[vert], cc)])^(-1 + length(all_neighbors(g, vert)))
    end

    d = Dict(edges(g) .=> weights)
    for e in edges(g)
        b = Int64(integrate(iv[(col[src(e)], col[dst(e)])] * cc))
        # b == 0 && continue
        for alph in 1:(b*d[e]-1)
            ans *= (alph * v.__attrs[:ev_dict][(col[src(e)], cc)] + (b * d[e] - alph) * v.__attrs[:ev_dict][(col[dst(e)], cc)]) // (b * d[e])
        end
    end

    return ans
end

function _push_omega(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, cc::CohomologyClass)::Cycle

    inters = integrate(beta * cc)
    try
        if inters < 2
            error("It is not a vector bundle")
        end
    catch e
        printstyled(stderr, "ERROR: ", bold=true, color=:red)
        printstyled(stderr, sprint(showerror, e), color=:light_red)
        println(stderr)
        return error_cycle()
    end

    return Cycle(Int64(inters - 1), 0)
end
