export
    R1_ev


"""
    R1_ev(l)

The equivariant class of the first derived functor of the pull-back of the toric line bundle ``l``.
# Arguments
- `l::ToricLineBundle`: the line bundle.

# Example
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,0}(\\mathbb{P}^{1},d)}\\mathrm{c_{top}}(R^{1}\\delta_{*}(\\mathrm{ev}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(-1)))^2 &= \\frac{1}{d^3} \\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using Oscar, ToricAtiyahBott)
julia> P1 = projective_space(NormalToricVariety, 1);

julia> beta = moment_graph(P1)[1,2];
(1,2) -> 1

julia> P = R1_ev(toric_line_bundle(P1, [-1]))^2;

julia> IntegrateAB(P1, beta, 0, P, show_bar=false);
Result: 1
```
"""
function R1_ev(l::ToricLineBundle)::EquivariantClass

    cc = cohomology_class(toric_divisor(l))

    rule = :(_R1_ev(v, od, nc, iv, g, col, weights, marks, $cc))

    return EquivariantClass(rule, eval(:((v, od, nc, iv, g, col, weights, marks) -> $rule)))
end

function _R1_ev(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Marks_type, cc::CohomologyClass)::T

    ans = F(1)

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
        for alph in (b*d[e]+1):-1
            ans *= (alph * v.__attrs[:ev_dict][(col[src(e)], cc)] + (b * d[e] - alph) * v.__attrs[:ev_dict][(col[dst(e)], cc)]) // (b * d[e])
        end
    end

    return ans
end

function _R1_ev(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, cc::CohomologyClass)::Cycle

    inters = integrate(beta * cc)
    try
        if inters > -1
            error("The R1 is not a vector bundle")
        end
    catch e
        printstyled(stderr, "ERROR: ", bold=true, color=:red)
        printstyled(stderr, sprint(showerror, e), color=:light_red)
        println(stderr)
        return error_cycle()
    end

    return Cycle(Int64(-inters - 1), 0)
end