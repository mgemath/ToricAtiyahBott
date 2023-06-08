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
```julia-repl
julia> v = projective_space(NormalToricVariety, 1);  # 1-dimensional projective space
julia> p = a_point(v); # the cohomology class of a point. Note that p^0 gives the class of the entire variety
julia> P = ev(1, p)*ev(2, p);
julia> IntegrateAB(v, p^0, 2, P);
Result: 1
julia> v = projective_space(NormalToricVariety, 2);  # 2-dimensional projective space
julia> l = toric_line_bundle(v, [ZZRingElem(1)]); 
julia> P = (ev(1, l)*ev(2, l))^2;
julia> line = cohomology_class(toric_divisor(v, [1,0,0]));
julia> IntegrateAB(v, line, 2, P);
Result: 1
```
!!! warning "Attention!"

    The program will stop if `j` is not between 1 and the number of marks.

"""
function ev(j::Int64, cc::Union{CohomologyClass,ToricLineBundle})::EquivariantClass
    
    SPLIT = split_cc(isa(cc, CohomologyClass) ? cc : cohomology_class(toric_divisor(cc)))
    Z = [(i, Int64.(exponents(SPLIT[i])[1,:]), Int64(Oscar.coefficients(SPLIT[i])[1])) for i in eachindex(SPLIT)]
    rule = :(_ev(v, od, nc, iv, g, col, weights, marks, $j, $Z))
    
    return EquivariantClass( rule, eval( :((v, od, nc, iv, g, col, weights, marks) -> $rule )))
end
# function ev(j::Int64, cc::CohomologyClass)::EquivariantClass
    
#     # rule = :(ev(v, od, nc, col, marks, $j, $Z))
#     # return EquivariantClass( rule, eval( :((v, od, nc, d, g, col, w, marks) -> $rule )))
#     rule = :(_ev(v, od, nc, iv, g, col, weights, marks, $j, $cc))
#     return EquivariantClass( rule, eval( :((v, od, nc, iv, g, col, weights, marks) -> $rule )))
# end

# function _ev(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, j::Int64, Z::CohomologyClass)::T

#     length(marks) == 0 && return one(od[first(keys(od))])
#     SPLIT = split_cc(Z)

#     ans = [one(od[first(keys(od))]) for _ in 1:length(SPLIT)]

#     for (i, mon) in enumerate(SPLIT)
#         e = exponents(mon)
#         for (k, ray) in enumerate(rays(v))
#             ####
#             if !(ray in rays(maximal_cones(v)[col[marks[j]]]))
#                 if e[1, k] == 0
#                     continue
#                 else
#                     # ans[i] *= zero(T)
#                     ans[i] *= zero(od[first(keys(od))])
#                     break
#                 end
#             end
#             ####
#             # e[1, k] == 0 && continue
#             # ray in rays(maximal_cones(v)[col[marks[j]]]) || continue
#             # l_j = e[1, k]
#             for n_gamma in nc[col[marks[j]]]
#                 ray in rays(maximal_cones(v)[n_gamma]) && continue
#                 ans[i] *= od[(col[marks[j]], n_gamma)]^Int64(e[1, k])
#                 break
#             end
#         end
#         ans[i] *= Int64(Oscar.coefficients(mon)[1])
#     end

#     return sum(ans)
# end


function _ev(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, j::Int64, Z::Vector{Tuple{Int64, Matrix{Int64}, Int64}})::T

    length(marks) == 0 && return one(od[first(keys(od))])
    ans = [one(od[first(keys(od))]) for _ in 1:length(Z)]

    for (i, e, coef) in Z
        for (k, ray) in enumerate(rays(v))
            ####
            if !(ray in rays(maximal_cones(v)[col[marks[j]]]))
                if e[k] == 0
                    continue
                else
                    # ans[i] *= zero(T)
                    ans[i] *= zero(od[first(keys(od))])
                    break
                end
            end
            ####
            # e[1, k] == 0 && continue
            # ray in rays(maximal_cones(v)[col[marks[j]]]) || continue
            # l_j = e[1, k]
            for n_gamma in nc[col[marks[j]]]
                ray in rays(maximal_cones(v)[n_gamma]) && continue
                ans[i] *= od[(col[marks[j]], n_gamma)]^e[k]
                break
            end
        end
        ans[i] *= coef
    end

    return sum(ans)
end

function _ev(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, i::Int64, Z::Vector{Tuple{Int64, Matrix{Int64}, Int64}})::Cycle

    try
        if !all(i -> sum(Z[1][2]) == sum(Z[i][2]), eachindex(Z)) #!ishomogeneous(polynomial(Z))
            error("The cohomology class is not homogeneous")
        end
        if (i < 1 || i > n_marks)
            error(string("ev requires a positive integer between 1 and ", n_marks, ", correct ",i))
        end
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    
    return Cycle(sum(Z[1][2]), 0)
end

# function ev_dict(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, j::Int64, Z::Vector{Tuple{Int64, Matrix{Int64}, Int64}})
#     ans = Dict{Int64, T}()

#     for c in 1:n_maximal_cones(v)
#         ans[c] = _ev(v, od, nc, iv, g, (c, ), weights, (1, ), 1, Z)
#     end

#     return ans
# end
# function __ev(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, j::Int64, Z::Vector{Tuple{Int64, Matrix{Int64}, Int64}}, evd)::T
#     return evd[col[marks[j]]]
# end