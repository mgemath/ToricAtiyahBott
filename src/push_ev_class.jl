export 
    push_ev

function push_ev(l::ToricLineBundle)::EquivariantClass

    cc = cohomology_class(toric_divisor(l))
    SPLIT = split_cc(cc)
    Z = [(i, Int64.(exponents(SPLIT[i])[1,:]), Int64(Oscar.coefficients(SPLIT[i])[1])) for i in eachindex(SPLIT)]

    rule = :(_push_ev(v, od, nc, iv, g, col, weights, marks, $cc, $Z))
    # rule = :(_push_ev(v, od, nc, iv, g, col, weights, marks, $l))
    
    return EquivariantClass( rule, eval( :((v, od, nc, iv, g, col, weights, marks) -> $rule )))
end

function _push_ev(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, cc::CohomologyClass, Z::Vector{Tuple{Int64, Matrix{Int64}, Int64}})::T

    ans = one(od[first(keys(od))])
    d = Dict(edges(g) .=> weights)

    evd = Dict{Int64, T}()
    for vert in 1:nv(g)
        if !(col[vert] in keys(evd))
            evd[col[vert]] = _ev(v, od, nc, iv, g, (col[vert], ), weights, (1, ), 1, Z)
        end

        1-length(all_neighbors(g, vert)) == 0 && continue
        evd[col[vert]] == 0 && return evd[col[vert]]
        ans *= (evd[col[vert]])^(1-length(all_neighbors(g, vert))) 
    end

    for e in edges(g)
        b = Int64( integrate(iv[(col[src(e)], col[dst(e)])]*cc) )
        b == 0 && continue
        for alph in 0:(b*d[e])
            ans *= (alph*evd[col[src(e)]]+(b*d[e]-alph)*evd[col[dst(e)]])//(b*d[e])
        end
    end

    return ans
end

function _push_ev(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, cc::CohomologyClass, Z::Vector{Tuple{Int64, Matrix{Int64}, Int64}})::Cycle

    inters = integrate(beta*cc)
    try
        if inters < 0
            error("The push-forward is not a vector bundle")
        end
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    
    return Cycle(Int64(inters + 1), 0)
end

# function _push_ev(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, l::ToricLineBundle)::T

#     ans = one(od[first(keys(od))])
#     d = Dict(edges(g) .=> weights)
#     SPLIT = split_cc( cohomology_class(toric_divisor(l)))
#     Z = [(i, Int64.(exponents(SPLIT[i])[1,:]), Int64(Oscar.coefficients(SPLIT[i])[1])) for i in eachindex(SPLIT)]

#     evd = Dict{Int64, T}()
#     for vert in 1:nv(g)
#         col[vert] in keys(evd) && continue
#         evd[col[vert]] = _ev(v, od, nc, iv, g, (col[vert], ), weights, (1, ), 1, Z)
#     end

#     for e in edges(g)
#         # sub_col = (col[src(e)], col[dst(e)])
#         # sub_mark = (1,2)
#         # e1 = _ev(v, od, nc, iv, g, sub_col, weights, sub_mark, 1, Z)
#         # e2 = _ev(v, od, nc, iv, g, sub_col, weights, sub_mark, 2, Z)

#         b = Int64( integrate(iv[(col[src(e)], col[dst(e)])]*cohomology_class(l)) )
#         b == 0 && continue
#         for alph in 0:(b*d[e])
#             ans *= (alph*evd[col[src(e)]]+(b*d[e]-alph)*evd[col[dst(e)]])//(b*d[e])
#         end
#         # for alph in 0:(b*d[e])
#         #     ans *= (alph*e1+(b*d[e]-alph)*e2)//(b*d[e])
#         # end
#         # for alph in 0:(b*d[e])
#         #     ans *= (alph*scal[col[src(e)]]+(b*d[e]-alph)*scal[col[dst(e)]])//d[e]
#         # end
#     end

#     for vert in 1:nv(g)
#         1-length(all_neighbors(g, vert)) == 0 && continue
#         # sub_col = (col[vert] ,)
#         # sub_mark = (1,)
#         # e1 = _ev(v, od, nc, iv, g, sub_col, weights, sub_mark, 1, Z)
#         # e1 == 0 && return e1

#         evd[col[vert]] == 0 && return evd[col[vert]]
#         # println("e1= ", e1)
#         # println("exp= ", 1-length(all_neighbors(g, vert)))
#         ans *= (evd[col[vert]])^(1-length(all_neighbors(g, vert)))   
#     end

#     return ans
# end