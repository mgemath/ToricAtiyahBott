export
    class_one


"""
    class_one()

Equivariant class equals to the number one. Useful to compute the degree of a moduli space of dimension 0.

# Example
Let ``v`` be a toric variety and ``\\beta`` be a class such that ``\\overline{M}_{0,0}(v,\\beta)`` is zero dimensional. The following Gromov-Witten invariants
```math
\\begin{equation}
\\int_{\\overline{M}_{0,0}(v,\\beta)}1 
\\end{equation}
```
can be computed, for example, as
```jldoctest; setup = :(using Oscar, ToricAtiyahBott)
julia> P3 = projective_space(NormalToricVariety, 3);

julia> x = domain(blow_up(P3, [1,1,1]; coordinate_name="Ex1"));

julia> v = domain(blow_up(x, [-1,0,0]; coordinate_name="Ex2"));

julia> mg = moment_graph(v, show_graph=false);

julia> (H, E1, E2) = (mg[7,8], mg[4,5], mg[1,2]);

julia> (d, e1, e2) = (2, -2, -2);

julia> beta = d*H + e1*E1 + e2*E2;

julia> P = class_one();

julia> IntegrateAB(v, beta, 0, P, show_bar=false);
Result: 1//8
```
"""
function class_one()::EquivariantClass

    rule = :(_class_one(v, od, nc, iv, g, col, weights, marks))

    return EquivariantClass(rule, eval(:((v, od, nc, iv, g, col, weights, marks) -> $rule)))
end

function _class_one(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}})::T

    return F(1)
end

function _class_one(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64)::Cycle

    return Cycle(0, 0)
end
