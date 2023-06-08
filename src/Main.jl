export 
    IntegrateAB

# T = QQFieldElem
# T = Rational{BigInt}
T = fpFieldElem
# T = AbstractAlgebra.Generic.Frac{QQMPolyRingElem}
# scal = [0]

"""
    IntegrateAB(v, beta, m, P; do_check)

Apply the Atiyah-Bott residue formula to compute the integral of the equivariant class `P` in the moduli space of rational marked stable maps to the toric variety `v` of class `beta` with `m` marks.
# Arguments
- `v::NormalToricVariety`: the toric variety.
- `beta::CohomologyClass`: the class of the stable maps, it must be the effective class of a curve.
- `m::Int64`: the number of marks.
- `P`: the equivariant class.
- `do_check::Bool`: if `true`, checks if `P` is a well defined zero cycle, and stops the computation if this is not true. If `false`, the computation may have an unexpected behaviour. By default is `true`.

In order to use this function, one must define `v`, `beta` and `P` using Oscar:
```julia-repl
julia> using Oscar
julia> v = del_pezzo_surface(1); # this is the blow-up of the projective plane at a point
julia> beta = cohomology_class(toric_divisor(v, [0,0,1,0])); # class of pull back of a line of P2
julia> P = ev(1, a_point(v))*ev(2, a_point(v)); # pull back of a point through the first and second evaluations maps
julia> IntegrateAB(v, beta, 2, P);
Result: 1//1
```

The function returns an array of the same dimension of `P` (non-vectorized classes are assumed as 1-dimensional arrays). The Julia notation for accessing to array is `name_of_array[i]` where `i` is an index starting from 1.

More examples are available in the support of the equivariant classes. It is enough to type `?` and then the name of the class. Currently, the supported classes are:

* `ev(j, cc)` (Euler class of ``\\mathrm{ev}_j^*(cc)``) where cc is either a cohomology class or a line bundle
* `Psi(a)`    (cycle of ``\\psi``-classes)
* `Jet(p, l)` (Euler class of the jet bundle ``J^p`` of the line bundle ``l``)

To add more classes, please contact the authors.
"""
function IntegrateAB(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, P_input; do_check::Bool = true)::Union{Vector{T},Nothing}
    
    if !issmooth(v)
        printstyled("ERROR: ", bold=true, color=:red)
        println("The variety must be smooth")
        return nothing
    end

    if !isprojective(v)
        printstyled("ERROR: ", bold=true, color=:red)
        println("The variety must be projective")
        return nothing
    end

    if !ishomogeneous(polynomial(beta))
        printstyled("ERROR: ", bold=true, color=:red)
        println("The class of the curve is not homogeneous")
        return nothing
    end

    if degree(Int64, polynomial(beta)) != dim(v)-1
        printstyled("ERROR: ", bold=true, color=:red)
        println("The class is not the class of a curve")
        return nothing
    end

    NEF = nef_divisors(v)
    if !is_effective_class(NEF, beta)
        printstyled("ERROR: ", bold=true, color=:red)
        println("The class of the curve is not effective")
        return nothing
    end

    # R, (x1, x2, x3) = polynomial_ring(QQ, ["x1", "x2", "x3"])
    # S = fraction_field(R)
    # (X1, X2, X3) = (S(x1), S(x2), S(x3))
    # s = ([X1, X2, X3]...,)

    # global T = QQFieldElem
    F = GF(2903)
    # F = QQ
    # s = (T.(rand(UInt16, nrays(v)))...,)

    s = (F.(rand(UInt16, nrays(v)))...,)

    # ans = zero(s[1])
    nc = neighbors_cones(v)
    od = omega_dict(v, nc, s)
    d = get_inv_curve(v, nc)
    otd::Dict{Int64, T} = Dict{Int64, T}() # omega_total dict
    for i in 1:n_maximal_cones(v)
        otd[i] = omega_total(od, nc, i)
    end
    
    max_n_vert = mapreduce(D -> Int64(integrate(D*beta)), +, NEF) + 1

    local n_results::Int64 = 1

    if isa(P_input, Array)
        n_results = length(P_input)
    end

    local P::Vector{Function} = Vector(undef, n_results)
    
    if isa(P_input, Array)
        for i in eachindex(P)
            P[i] = P_input[i].func
        end
    else
        P[1] = P_input.func
    end
    
    if do_check && !is_zero_cycle(v, beta, n_marks, P)
        return nothing
    end
    
    local result::Vector{Vector{T}} = [[zero(s[1]) for _ in 1:n_results] for _ in 1:Threads.nthreads()]
    local partial_res::Vector{T} = [zero(s[1]) for _ in 1:1:n_results]

    # P = P_input.func

    for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert])
        g = LStoGraph(ls)
        # top_aut::Int64 = count_iso(ls)
        (rev_dfs, parents, has_ci, subgraph_ends_rev, subgraph_ends, left_siblings) = init_ColorsIt(ls)
        for col in ColorsIt(ls, n_maximal_cones(v), rev_dfs, parents, has_ci, subgraph_ends_rev, subgraph_ends, left_siblings)
            c = Int64.(col)
            is_admissible_color(nc, g, c) || continue
            top_aut::Int64 = count_iso(ls, col)
            
            for w in multip(Ms(v, NEF, g, c, d, beta))

                local store_aut::Dict{Vector{Int64},Int64} = Dict{Vector{Int64},Int64}()
                local temp_m::Vector{Int64}
                local aut::Int64

                #for m in Base.Iterators.product(repeat([1:nv(g)], n_marks)...)    #we run among all marks of g, if n_marks==0 we have only the empty mark                           
                for m in Base.Iterators.filter(m -> top_aut == 1 || isempty(m) || maximum(m) < 3 || ismin(ls, col, m, parents, subgraph_ends), Base.Iterators.product(repeat([1:nv(g)], n_marks)...))
                    for res in eachindex(partial_res)
                        partial_res[res] = Base.invokelatest(P[res], v, od, nc, d, g, c, w, m)
                    end
                    # foreach(res -> eq!(partial_res[res], Base.invokelatest(P[res], g,c,w,s,m)), eachindex(partial_res))
                    all(res -> partial_res[res] == zero(s[1]), eachindex(partial_res)) && continue # check if at least
                    # eq == 0 && continue
                    
                    # col_rel(v,w) = (v == w ||(c[v] == c[w] && !(v in m) && !(w in m)))
                    # local aut::Int64 = Graphs.Experimental.count_isomorph(g,g,vertex_relation=col_rel)
                    temp_m = sort(unique(m))
                    
                    if !haskey(store_aut, temp_m)
                        store_aut[temp_m] = count_iso(ls, col, m)
                    end
                    aut = store_aut[temp_m]


                    E = Euler_inv(v, od, nc, otd, d, g, c, w, m)//(aut*prod(w))
                    for res in eachindex(partial_res)      # compute each term of the array P
                        partial_res[res] *= E
                        # println(partial_res[res])
                        result[Threads.threadid()][res] += partial_res[res]
                    end
                    # for w in multip(Ms(v, NEF, g, c, d, beta))
                        #P = ev_j(X, od, nc, c, m, 1, point)*ev_j(X, od, nc, c, m, 2, point) #X, r_e, s, 1, c, m, point)*ev_j(X, r_e, s, 2, c, m, point)*ev_j(X, r_e, s, 3, c, m, point)*ev_j(X, r_e, s, 4, c, m, point)*ev_j(X, r_e, s, 5, c, m, point)
                        # P = ev_j(X, r_e, s, 1, c, m, line)*ev_j(X, r_e, s, 2, c, m, line)*ev_j(X, r_e, s, 3, c, m, line)*ev_j(X, r_e, s, 4, c, m, line)*ev_j(X, r_e, s, 5, c, m, line)
                    # ans += eq*E
                end
            end
        end
    end

    for nt in 2:Threads.nthreads()
        result[1] += result[nt]
    end
    
    if n_results == 1
        println("Result: ", result[1][1])
    else 
        for res in 1:n_results
            println("Result number ", res, ": ", result[1][res])
        end
    end
    return result[1]
end