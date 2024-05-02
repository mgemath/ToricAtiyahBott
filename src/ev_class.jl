export
    ev

"""
    ev(j, cc)
    ev(j, l)

Equivariant class of the pull-back of the cohomology class ``cc`` (or the toric line bundle ``l``) with respect to the j-th evaluation map.
# Arguments
- `j::Int64`: the evaluation map.
- `cc::CohomologyClass`: the cohomology class.
- `l::ToricLineBundle`: the line bundle.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{1},1)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{1}}(1)\\cdot\\mathrm{ev}_{2}^{*}\\mathcal{O}_{\\mathbb{P}^{1}}(1) &= 1 \\\\
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{2},1)}(\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)\\cdot\\mathrm{ev}_{2}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1))^2 &= 1 \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using Oscar, ToricAtiyahBott)
julia> v = projective_space(NormalToricVariety, 1);  # 1-dimensional projective space

julia> p = a_point(v); # the cohomology class of a point. Note that p^0 gives the class of the entire variety

julia> P = ev(1, p)*ev(2, p);

julia> IntegrateAB(v, p^0, 2, P, show_bar=false); # show_bar can be also true
Result: 1

julia> v = projective_space(NormalToricVariety, 2);  # 2-dimensional projective space

julia> l = toric_line_bundle(v, [1]);

julia> P = (ev(1, l)*ev(2, l))^2;

julia> line = cohomology_class(toric_divisor(v, [1,0,0]));

julia> IntegrateAB(v, line, 2, P, show_bar=false); # show_bar can be also true
Result: 1
```
!!! warning "Attention!"

    The program will stop if `j` is not between 1 and the number of marks.


Let us give some more examples. Let ``v = \\mathbb{P}(\\mathcal{O}_{\\mathbb{P}^3}\\oplus\\mathcal{O}_{\\mathbb{P}^3}(5))``.
```julia-repl
julia> P3 = projective_space(NormalToricVariety, 3);
julia> v = proj(toric_line_bundle(P3, [0]),toric_line_bundle(P3, [5]));
```
Using [`moment_graph`](@ref) we have a quick access to the moment graph of ``v``.
```julia-repl
julia> mg = moment_graph(v);
```
Consider the following curve class.
```julia-repl
julia> beta = mg[1,2];
```
If ``\\mathrm{p}`` is the class of a point of ``v``, in order to compute the following invariant
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,1}(v, \\beta)} \\mathrm{ev}_{1}^{*}(\\mathrm{p}) &= 1 \\\\
\\end{aligned}
```
we use the code:
```julia-repl
julia> P = ev(1, a_point(v));
julia> IntegrateAB(v, beta, 1, P);
```
In order to speed up the computation, many equivariant classes of the same moduli space can be vectorized. For example the following two invariants are in the same moduli space.
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^2, 1)} \\mathrm{ev}_{1}^{*}(\\mathrm{p})\\cdot\\mathrm{ev}_{2}^{*}(\\mathrm{p}) &= 1 \\\\
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^2, 1)} \\mathrm{ev}_{1}^{*}(\\mathrm{p})\\cdot\\psi_{1}^{1}\\cdot\\psi_{2}^{1} &= -1. \\
\\end{aligned}
```
The best way to compute them is by defining an array `P=[P1,P2]` and compute them together.
```jldoctest; setup = :(using Oscar, ToricAtiyahBott)
julia> v = projective_space(NormalToricVariety, 2);

julia> p = a_point(v);

julia> P1 = ev(1, p)*ev(2, p);

julia> P2 = ev(1, p)*Psi([1,1]);

julia> P = [P1,P2];

julia> line = cohomology_class(toric_divisor(v, [1,0,0]));

julia> IntegrateAB(v, line, 2, P, show_bar=false);
Result number 1: 1
Result number 2: -1
```
"""
function ev(j::Int64, cc::CohomologyClass)::EquivariantClass

    rule = :(_ev(v, od, nc, iv, g, col, weights, marks, $j, $cc))
    return EquivariantClass(rule, eval(:((v, od, nc, iv, g, col, weights, marks) -> $rule)))
end

function ev(j::Int64, l::ToricLineBundle)::EquivariantClass

    return ev(j, cohomology_class(l))
end

function _ev(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Marks_type, j::Int64, cc::CohomologyClass)::T

    length(marks) == 0 && return F(1)
    if (col[marks[j]], cc) in keys(v.__attrs[:ev_dict])
        return v.__attrs[:ev_dict][(col[marks[j]], cc)]
    end

    v.__attrs[:ev_dict][(col[marks[j]], cc)] = just_ev(v, od, nc, col[marks[j]], cc)

    return v.__attrs[:ev_dict][(col[marks[j]], cc)]
end

function just_ev(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, col_vert::Int64, cc::CohomologyClass)::T

    SPLIT = split_cc(cc)
    Z = [(i, Int64.(exponents(SPLIT[i])[1, :]), Oscar.coefficients(SPLIT[i])[1]) for i in eachindex(SPLIT)]
    ans = [F(1) for _ in 1:length(Z)]

    for (i, e, coef) in Z
        for (k, ray) in enumerate(rays(v))
            e[k] == 0 && continue

            if !(ray in rays(maximal_cones(v)[col_vert]))
                ans[i] = F(0)
                break
            end

            for n_gamma in nc[col_vert]
                ray in rays(maximal_cones(v)[n_gamma]) && continue
                # ans[i] *= od[(col_vert, n_gamma)]^e[k]
                mul!(ans[i], ans[i], od[(col_vert, n_gamma)]^e[k])
                break
            end
        end

        ans[i] *= coef
    end

    return sum(ans)
end


function _ev(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, j::Int64, cc::CohomologyClass)::Cycle

    SPLIT = split_cc(cc)
    Z = [(i, Int64.(exponents(SPLIT[i])[1, :]), Oscar.coefficients(SPLIT[i])[1]) for i in eachindex(SPLIT)]

    try
        if !all(i -> sum(Z[1][2]) == sum(Z[i][2]), eachindex(Z)) #!ishomogeneous(polynomial(Z))
            error("The cohomology class is not homogeneous")
        end
        if (j < 1 || j > n_marks)
            error(string("ev requires a positive integer between 1 and ", n_marks, ", correct ", j))
        end
    catch e
        printstyled(stderr, "ERROR: ", bold=true, color=:red)
        printstyled(stderr, sprint(showerror, e), color=:light_red)
        println(stderr)
        return error_cycle()
    end

    return Cycle(sum(Z[1][2]), 0)
end
