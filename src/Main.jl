export
    IntegrateAB


"""
    IntegrateAB(v, beta, m, P; do_check, show_bar)

Apply the Atiyah-Bott residue formula to compute the integral of the equivariant class `P` in the moduli space of rational marked stable maps to the toric variety `v` of class `beta` with `m` marks.
# Arguments
- `v::NormalToricVariety`: the toric variety.
- `beta::CohomologyClass`: the class of the stable maps, it must be the effective class of a curve.
- `m::Int64`: the number of marks.
- `P`: the equivariant class.
- `do_check::Bool`: if `true`, checks if `P` is a well defined zero cycle, and stops the computation if this is not true. If `false`, the computation may have an unexpected behaviour. By default is `true`.
- `show_bar::Bool`: hide the progress bar if and only if this condition is `false`. By default is `true`.

In order to use this function, one must define `v`, `beta` and `P` using Oscar:
```jldoctest; setup = :(using Oscar, ToricAtiyahBott)
julia> v = del_pezzo_surface(NormalToricVariety, 1); # this is the blow-up of the projective plane at a point

julia> beta = cohomology_class(toric_divisor(v, [0,0,1,0])); # class of pull back of a line of P2

julia> P = ev(1, a_point(v))*ev(2, a_point(v)); # pull back of a point through the first and second evaluations maps

julia> IntegrateAB(v, beta, 2, P, show_bar=false); # show_bar can be also true
Result: 1
```
The function returns an array of the same dimension of `P` (non-vectorized classes are assumed as 1-dimensional arrays). The Julia notation for accessing to array is `name_of_array[i]` where `i` is an index starting from 1.

More examples are available in the support of the equivariant classes. It is enough to type `?` and then the name of the class. Currently, the supported classes are:

* `ev(j, cc)` (Euler class of ``\\mathrm{ev}_j^*(cc)`` where ``cc`` is either a cohomology class or a line bundle)
* `Psi(a)`    (cycle of ``\\psi``-classes)
* `push_ev(l)`(push-forward under the forget map of ``\\mathrm{ev}_j^*(l)``)
* `R1_ev(l)`  (first derived functor of the pull-back of ``\\mathrm{ev}_j^*(l)``)
* `Jet(p, l)` (Euler class of the jet bundle ``J^p`` of the line bundle ``l``)

To add more classes, please contact the authors.
"""
function IntegrateAB(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, P_input; do_check::Bool=true, show_bar::Bool=true)::Union{Vector{T},Nothing}

    if !is_smooth(v)
        printstyled("ERROR: ", bold=true, color=:red)
        println("The variety must be smooth")
        return nothing
    end

    if !is_projective(v)
        printstyled("ERROR: ", bold=true, color=:red)
        println("The variety must be projective")
        return nothing
    end

    if !is_homogeneous(polynomial(beta))
        printstyled("ERROR: ", bold=true, color=:red)
        println("The class of the curve is not homogeneous")
        return nothing
    end

    if degree(Int64, polynomial(beta)) != dim(v) - 1
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

    s = (F.(rand(UInt16, n_rays(v)))...,)

    nc = neighbors_cones(v)
    od = omega_dict(v, nc, s)
    d = get_inv_curve(v, nc)
    otd::Dict{Int64,T} = Dict{Int64,T}() # omega_total dict
    for i in 1:n_maximal_cones(v)
        otd[i] = omega_total(od, nc, i)
    end

    v.__attrs[:ev_dict] = Dict{Tuple{Int64,CohomologyClass},T}()

    local n_results::Int64 = 1

    if isa(P_input, Array)
        n_results = length(P_input)
    end

    local P::Vector{Function} = Vector{Function}(undef, n_results)

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

    local result::Vector{T} = [F(0) for _ in 1:n_results]
    local partial_res::Vector{T} = [F(0) for _ in 1:n_results]

    max_n_vert = mapreduce(D -> Int64(integrate(D * beta)), +, NEF, init=1)
    if show_bar #set up progress data
        number_trees = A000055(max_n_vert)
        # threshold = sum(vert -> number_trees[vert] * (n_maximal_cones(v)) * ((length(nc[1]))^(vert - 1)) * (vert^n_marks), 2:max_n_vert)
        threshold = sum(vert -> number_trees[vert] * (n_maximal_cones(v)) * ((length(nc[1]))^(vert - 1)) * binomial(vert+n_marks-1, n_marks), 2:max_n_vert)
        progress_bar::Progress = Progress(threshold, barglyphs=BarGlyphs("[=> ]"), color=:green)
        current_graph::Threads.Atomic{Int64} = Threads.Atomic{Int}(0)
    end

    Lambda_gamma_e_dict::Dict{Tuple{Int64, Int64, Int64}, T} = Dict{Tuple{Int64, Int64, Int64}, T}()
    
    for ls in Iterators.flatten([TreeIt(i) for i in 2:max_n_vert])

        tree_aut = count_iso(ls)
        g = LStoGraph(ls)

        CI, parents, subgraph_ends = col_it_init(ls, nc)
        for col in CI
            top_aut::Int64 = count_iso(ls, col)

            aut::Int64 = 0

            MULTI = collect(multip(Ms(v, NEF, g, col, d, beta)))
            # last_seen_unique = Int64[]

            # for m_inv in Base.Iterators.filter(m -> top_aut == 1 || isempty(m) || maximum(m) < 3 || ismin(ls, col, (m...,), parents, subgraph_ends), with_replacement_combinations(1:nv(g), n_marks))
            for m_inv in with_replacement_combinations(1:nv(g), n_marks)
                aut = count_iso(ls, col, m_inv)
                
                # u = unique(m_inv)
                # if last_seen_unique != u
                #     last_seen_unique = u
                #     aut = count_iso(ls, col, (m_inv...,))
                # end
                
                for w in MULTI
                    PRODW = prod(w)
                    E = F(0)
                    D = Dict(edges(g) .=> w)
                    # for m_v in multiset_permutations(m_inv, n_marks)
                    for m in Base.Iterators.filter(mul_per -> top_aut == 1 || isempty(mul_per) || maximum(mul_per) < 3 || ismin(ls, col, mul_per, parents, subgraph_ends), multiset_permutations(m_inv, n_marks))

                        for res in eachindex(partial_res)
                            partial_res[res] = Base.invokelatest(P[res], v, od, nc, d, g, col, w, m)
                        end

                        all(res -> partial_res[res] == F(0), eachindex(partial_res)) && continue # check if at least one partial result is not zero
                   
                        if E == F(0)
                            E = Euler_inv(v, od, nc, otd, d, g, col, w, m) // (aut * PRODW)
                            for e in edges(g)
                                triple = (D[e], min(col[src(e)], col[dst(e)]), max(col[src(e)], col[dst(e)]))
                                if !haskey(Lambda_gamma_e_dict, triple)
                                    Lambda_gamma_e_dict[triple] = Lambda_Gamma_e(v, od, nc, d, triple[1], triple[2], triple[3])
                                end
                                mul!(E, E, Lambda_gamma_e_dict[triple])
                            end
                        end

                        for res in eachindex(partial_res)      # compute each term of the array P
                            partial_res[res] *= E
                            result[res] += partial_res[res]
                        end
                    end
                end
                if show_bar #update the progress bar
                    # Threads.atomic_add!(current_graph, tree_aut ÷ aut)
                    Threads.atomic_add!(current_graph, tree_aut ÷ top_aut)
                    # Threads.lock(l)
                    update!(progress_bar, current_graph[],
                        showvalues=[(:"Total number of graphs", threshold), (:"Current graph", current_graph[])])
                    # Threads.unlock(l)
                end
            end
        end

    end

    if n_results == 1
        println("Result: ", result[1])
    else
        for res in 1:n_results
            println("Result number ", res, ": ", result[res])
        end
    end
    return result[:]
end