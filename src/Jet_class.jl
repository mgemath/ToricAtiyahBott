export
    Jet


"""
    Jet(p, l)

Equivariant class of the jet bundle ``J^p`` of the pull back of the toric line bundle ``l`` with respect to the first ``\\psi``-class.
# Arguments
- `p::Int64`: the exponent of the Jet bundle. In particular, it is a bundle of rank ``p+1``.
- `l::Int64`: the toric line bundle that is pulled back.


!!! note

    In order to define this bundle, the number of marks must be at least 1.
    You cannot multiply this bundle by the class `Psi(a)`.


# Example
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{3},d)}\\frac{\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(1)^{2}}{k}\\cdot\\mathrm{c_{top}}(J^{4d-2}(\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{3}}(k))) &= \\frac{(4d-2)!}{(d!)^{4}} \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{2},1)}\\mathrm{c_{top}}(J^{2}(\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(3))) &= 9 \\\\
\\end{aligned}
```
can be computed as
```jldoctest; setup = :(using Oscar, ToricAtiyahBott)
julia> d=1;k=1; #for other values of d, change this line

julia> v = projective_space(NormalToricVariety, 3);

julia> line = cohomology_class(toric_divisor(v, [1,0,0,0]))^2;

julia> beta = d*line;

julia> l = toric_line_bundle(v, [k]);

julia> P = ev(1,line)//k*Jet(4*d-2,l);

julia> IntegrateAB(v, beta, 1, P, show_bar=false);   #The value of this integral does not depend on k, only on d
Result: 2

julia> v = projective_space(NormalToricVariety, 2);

julia> l = toric_line_bundle(v, [1]);

julia> line = cohomology_class(l);

julia> IntegrateAB(v, line, 1, Jet(2, l^3), show_bar=false); # flexes of a cubic curve
Result: 9
```
"""
function Jet(p::Int64, l::ToricLineBundle)::EquivariantClass

    stirling = stirling_tuple(p + 1)
    cc = cohomology_class(toric_divisor(l))

    rule = :(_Jet(v, od, nc, iv, g, col, weights, marks, $p, $cc, $stirling))

    return EquivariantClass(rule, eval(:((v, od, nc, iv, g, col, weights, marks) -> $rule)))
end

function _Jet(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, p::Int64, cc::CohomologyClass, st::Tuple{Vararg{Int64}})::Cycle

    try
        1 > n_marks && error("Jet is defined when there is at least one mark")
        0 > p && error(string("exponents of Jet must be nonnegative, correct ", p))
    catch e
        printstyled(stderr, "ERROR: ", bold=true, color=:red)
        printstyled(stderr, sprint(showerror, e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    psi_deg = p == 0 ? 0 : 1
    return Cycle(p + 1, psi_deg)
end

function _Jet(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Marks_type, p::Int64, cc::CohomologyClass, st::Tuple{Vararg{Int64}})::T

    local ans = F(0)

    if !((col[marks[1]], cc) in keys(v.__attrs[:ev_dict]))
        v.__attrs[:ev_dict][(col[marks[1]], cc)] = just_ev(v, od, nc, col[marks[1]], cc)
    end

    e = v.__attrs[:ev_dict][(col[marks[1]], cc)]


    if e != 0
        # temp = F(1)
        for h in p:-1:0
            # mul!(temp, temp, e)
            # add!(ans, ans, st[p+1-h]*(e^(p+1-h))*_Psi(v, od, nc, iv, g, col, weights, marks, [h]))
            ans += st[p+1-h] * (e^(p + 1 - h)) * _Psi(v, od, nc, iv, g, col, weights, marks, [h])
        end
    end

    return ans
end
