export
    push_ev


"""
    push_ev(l)

Equivariant class of the push-forward under the forget map of the pull-back of the  toric line bundle ``l``.
# Arguments
- `l::ToricLineBundle`: the line bundle.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{4},1)}\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{4}}(5))) &= 2875 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{4},1)}\\mathrm{c_{top}}(\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{4}}(5)))\\cdot\\psi_{1}^{1} &= -5750
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using Oscar, ToricAtiyahBott)
julia> v = projective_space(NormalToricVariety, 4);

julia> line = cohomology_class(toric_line_bundle(v, [1]))^3;

julia> P = push_ev(toric_line_bundle(v, [5]));

julia> IntegrateAB(v, line, 0, P, show_bar=false);
Result: 2875

julia> P = push_ev(toric_line_bundle(v, [5]))*Psi(1);

julia> IntegrateAB(v, line, 1, P, show_bar=false);
Result: -5750
```
"""
function push_ev(l::ToricLineBundle)::EquivariantClass

    cc = cohomology_class(toric_divisor(l))
    # SPLIT = split_cc(cc)
    # Z = [(i, Int64.(exponents(SPLIT[i])[1,:]), Int64(Oscar.coefficients(SPLIT[i])[1])) for i in eachindex(SPLIT)]

    rule = :(_push_ev(v, od, nc, iv, g, col, weights, marks, $cc))
    # rule = :(_push_ev(v, od, nc, iv, g, col, weights, marks, $l))

    return EquivariantClass(rule, eval(:((v, od, nc, iv, g, col, weights, marks) -> $rule)))
end

# function _push_ev(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, cc::CohomologyClass, Z::Vector{Tuple{Int64, Matrix{Int64}, Int64}})::T
function _push_ev(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, cc::CohomologyClass)::T
    ans = F(1)

    # evd = Dict{Int64, T}()

    for vert in 1:nv(g)
        # if !(col[vert] in keys(evd))
        #     evd[col[vert]] = _ev(v, od, nc, iv, g, (col[vert], ), weights, (1, ), 1, Z)
        # end
        if !((col[vert], cc) in keys(v.__attrs[:ev_dict]))
            v.__attrs[:ev_dict][(col[vert], cc)] = just_ev(v, od, nc, col[vert], cc)
        end

        1 - length(all_neighbors(g, vert)) == 0 && continue
        v.__attrs[:ev_dict][(col[vert], cc)] == 0 && return v.__attrs[:ev_dict][(col[vert], cc)]
        ans *= (v.__attrs[:ev_dict][(col[vert], cc)])^(1 - length(all_neighbors(g, vert)))
    end

    d = Dict(edges(g) .=> weights)
    for e in edges(g)
        b = Int64(integrate(iv[(col[src(e)], col[dst(e)])] * cc))
        b == 0 && continue
        for alph in 0:(b*d[e])
            ans *= (alph * v.__attrs[:ev_dict][(col[src(e)], cc)] + (b * d[e] - alph) * v.__attrs[:ev_dict][(col[dst(e)], cc)]) // (b * d[e])
        end
    end

    return ans
end

function _push_ev(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, cc::CohomologyClass)::Cycle

    inters = integrate(beta * cc)
    try
        if inters < 0
            error("The push-forward is not a vector bundle")
        end
    catch e
        printstyled(stderr, "ERROR: ", bold=true, color=:red)
        printstyled(stderr, sprint(showerror, e), color=:light_red)
        println(stderr)
        return error_cycle()
    end

    return Cycle(Int64(inters + 1), 0)
end
