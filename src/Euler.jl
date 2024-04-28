function Euler_inv(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, otd::Dict{Int64,T}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}})::T

    d = Dict(edges(g) .=> weights)

    ans = F(1)

    for e in edges(g)
        # ans *= Lambda_Gamma_e(v, od, nc, iv, col, d[e], e)
        mul!(ans, ans, Lambda_Gamma_e(v, od, nc, iv, col, d[e], e))
    end

    for O in 1:nv(g)
        t_O = length(neighbors(g, O))  # number of edges at O
        s_O = num_marks(marks, O)     # number of marks at O

        if t_O != 1
            # ans *= otd[col[O]]^(t_O - 1)
            mul!(ans, ans, otd[col[O]]^(t_O-1))
        end

        SUM = F(1)

        if t_O + s_O - 3 != 0
            SUM = F(0)
            for W in neighbors(g, O)
                e = Edge(max(O, W), min(O, W))
                # SUM += (omega_F(od, col[O], col[W], d[e])^(-1))
                add!(SUM, SUM, omega_F(od, col[O], col[W], d[e])^(-1))
            end
        end

        # PROD = one(T)
        PROD = F(1)
        for W in neighbors(g, O)
            e = Edge(max(O, W), min(O, W))#e = Edge(min(O,W),max(O,W))

            # PROD *= omega_F(od, col[O], col[W], d[e])
            mul!(PROD, PROD, omega_F(od, col[O], col[W], d[e]))

        end

        # ans *= (SUM^(t_O + s_O - 3)) // PROD
        mul!(ans, ans, (SUM^(t_O+s_O-3))//PROD)
    end

    return ans
end

function Lambda_Gamma_e(v::NormalToricVariety, od::Dict{Tuple{Int64,Int64},T}, nc::Dict{Int64,Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, col::Tuple{Vararg{Int64}}, d_e::Int64, e::Edge) #

    n_sigma2 = col[src(e)]
    n_sigma1 = col[dst(e)]

    b = od[(n_sigma1, n_sigma2)] // d_e #omega(v, n_sigma1, n_sigma2, scalars)//d_e

    ans = ((-1)^d_e) // ((factorial(d_e)^2) * (b^(2 * d_e)))
    # ans = (((-1)^d_e)*(d_e^(2*d_e))) // ((factorial(d_e)^2)*(omega(v, n_sigma1, n_sigma2, scalars)^(2*d_e)))

    # for gamma in maximal_cones(v) # all colors, different from sigma1 and sigma2, that meet sigma1
    for n_gamma in nc[n_sigma1]
        n_gamma == n_sigma2 && continue

        lambda_gamma_e = d_e

        # Now we have a gamma as above. Let us compute lambda_gamma_e
        # let us find the ray that is in sigma1 but not in gamma
        rays_gamma = rays(maximal_cones(v)[n_gamma])
        index = findfirst(ray -> !(ray in rays_gamma), rays(maximal_cones(v)[n_sigma1]))
        ray = rays(maximal_cones(v)[n_sigma1])[index]
        # for ray in rays(maximal_cones(v)[n_sigma1]) # let us find the ray that is in sigma1 but not in gamma
        #     ray in rays(maximal_cones(v)[n_gamma]) && continue

            p = [0 for _ in 1:n_rays(v)]
            p[findfirst(j -> j == ray, rays(v))] = 1

            class_e = iv[(n_sigma1, n_sigma2)] #[C_sigma1,sigma2]

            lambda_gamma_e *= Int64(integrate(class_e * cohomology_class(toric_divisor(v, p)))) # d_e * [C_sigma1,sigma2]*R_rj
        #     break
        # end

        lambda_gamma_e == -1 && continue

        i_range::UnitRange{Int64} = 0:lambda_gamma_e
        exp::Int64 = -1

        if lambda_gamma_e <= -2
            i_range = (lambda_gamma_e+1):-1
            exp = 1
        end

        a = od[(n_sigma1, n_gamma)]

        for i in i_range
            # ans *= (a - i * b)^exp
            mul!(ans, ans, (a-i*b)^exp)
        end
    end

    return ans
end