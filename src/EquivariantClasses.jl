export
    Psi,
    Jet

##################
###### Psi #######
##################

"""
    Psi(a)

Equivariant class of the cycle of ``\\psi``-classes.
# Arguments
- `a::Vector{Int64}`: the vector of the exponents of the ``\\psi`` classes. It is ordered, meaning that the first element is the exponent of ``\\psi_1``, the second is the exponent of ``\\psi_2``, and so on. Alternatively, `a` can be a number. In this case it is equivalent to [1,0,0,...,0].

!!! note

    The size of `a` must be at most `m`. If it is smaller, missing exponents will be considered as zeros.
    If `a` is a number, it will be considered as the exponent of ``\\psi_1``.

!!! warning "Attention!"

    The program will stop if we have one of the following conditions:

    * the size of `a` is bigger than `m`,
    * `a` contains a negative number.

# Example
The following Gromov-Witten invariants
```math
\\begin{aligned}
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{2},2)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2}\\cdot\\psi_{1}^{4} &= \\frac{1}{8} \\\\
\\int_{\\overline{M}_{0,2}(\\mathbb{P}^{2},1)}\\psi_{1}^{2}\\psi_{2}^{2} &= 6 \\\\
\\int_{\\overline{M}_{0,1}(\\mathbb{P}^{3},2)}\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)\\cdot(\\psi_{1}^{7}\\cdot\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)+\\psi_{1}^{6}\\cdot\\mathrm{ev}_{1}^{*}\\mathcal{O}_{\\mathbb{P}^{2}}(1)^{2}) &= -\\frac{5}{16} \\\\
\\end{aligned}
```
can be computed as
```julia-repl
julia> v = projective_space(NormalToricVariety, 2);  # 2-dimensional projective space
julia> line = cohomology_class(toric_divisor(v, [1,0,0])); # the cohomology class of a line
julia> P = ev(1, line)^2*Psi(4);
julia> IntegrateAB(v, 2*line, 1, P);
Result: 1//8
julia> Q = Psi([2,2])
julia> IntegrateAB(v, line, 2, Q);
Result: 6//1
julia> v = projective_space(NormalToricVariety, 3);  # 3-dimensional projective space
julia> plane = cohomology_class(toric_divisor(v, [1,0,0,0])); # the cohomology class of a plane
julia> P = ev(1, plane)*(Psi(7)*ev(1, plane)+Psi(6)*ev(1, plane)^2);
julia> IntegrateAB(v, 2*plane^2, 1, P);
Result: -5//16
```
!!! warning "Psi is singleton!"

    `Psi` cannot be multiplied by itself.
    ```julia-repl
    julia> v = projective_space(NormalToricVariety, 2);  # 2-dimensional projective space
    julia> line = cohomology_class(toric_divisor(v, [1,0,0])); # the cohomology class of a line
    julia> P = ev(1, line)^2*Psi(1)^4;                  #this is **wrong**
    julia> IntegrateAB(v, 2*line, 1, P);
    Warning: more instances of Psi has been found. Type:
    julia> ?Psi
    for support.
    julia> P = ev(1, line)^2*Psi(3)*Psi(1);            #this is **wrong**
    julia> IntegrateAB(v, 2*line, 1, P);
    Warning: more instances of Psi has been found. Type:
    julia> ?Psi
    for support.
    julia> P = ev(1, line)^2*Psi(4);
    julia> IntegrateAB(v, 2*line, 1, P);
    Result: 1//8
    ```
"""
function Psi(a)::EquivariantClass
    
    rule = :(_Psi(v, od, nc, iv, g, col, weights, marks, $a))
    return EquivariantClass( rule, eval( :((v, od, nc, iv, g, col, weights, marks) -> $rule )))
end

function _Psi(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, a::Int64)::T
    return _Psi(v, od, nc, iv, g, col, weights, marks, [a])
end


function _Psi(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, a::Vector{Int64})::T

    findfirst(x -> x>0, a) === nothing && return one(od[first(keys(od))]) #if all of them are zero or a is empty
    
    ans = one(od[first(keys(od))])

    # local q1::fmpq = fmpq(1)
    # local temp1::fmpq = fmpq(1)
    local Sum_ai::Int64
    local n::Int64
    local M::Int64
    local d = Dict(edges(g).=> weights) #assign weights to edges
    local inv_marks::Dict{Int64,Vector{Int64}} = invert_marks(marks, nv(g))
    
    for v in 1:nv(g)
        
        a_v = Int64[]
        for i in inv_marks[v]
            (i > length(a) || a[i] == 0) && continue
            push!(a_v, a[i])
        end

        Sum_ai = sum(a_v)
        Sum_ai == 0 && continue #if S contains only zeros, or it is empty, continue

        n = length(all_neighbors(g, v)) + length(inv_marks[v])
        
        n > 2 && Sum_ai > n - 3 && return zero(od[first(keys(od))])
        
        #If no previous condition holds, then n>1
        if n == 2 #necessary |S_v| == 1
            M = (-1)^a_v[1]
        else # n>2 and Sum_ai <= n - 3
            M = multinomial((n - 3 - Sum_ai, a_v...,))
        end
        

        local s1 = zero(od[first(keys(od))])
        
        for w in all_neighbors(g, v)
            s1 += d[Edge(max(v,w),min(v,w))]//od[col[v],col[w]]
        end
        ans *= M*(s1^(-Sum_ai))
    end
        
    return ans
end

function _Psi(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, a::Int64)::Cycle
    
    try
        a < 0 && error(string("exponents of psi classes must be nonnegative, correct ", a))
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    psi_deg = a==0 ? 0 : 1

    return Cycle(a, psi_deg)
end

function _Psi(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, a::Vector{Int64})::Cycle
    
    try
        length(a) > n_marks && error(string("size of ",a," is greater than ",n_marks))

        if findfirst(x -> x<0, a) !== nothing #if some of them is negative
            error(string("exponents of psi classes must be nonnegative, correct ", a))
        end
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    psi_deg = sum(a)==0 ? 0 : 1
    return Cycle(sum(a), psi_deg)
end

##################
###### Jet #######
##################
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
\\end{aligned}
```
can be computed as
```julia-repl
julia> d=1;k=1; #for other values of d, change this line
julia> v = projective_space(NormalToricVariety, 3);
julia> line = cohomology_class(toric_divisor(v, [1,0,0,0]))^2;
julia> beta = d*line;
julia> l = toric_line_bundle(v, [ZZRingElem(k)]);
julia> P = ev(1,line)//k*Jet(4*d-2,l);
julia> IntegrateAB(v, beta, 1, P);   #The value of this integral does not depend on k, only on d
```
"""
function Jet(p::Int64, l::ToricLineBundle)::EquivariantClass

    stirling = stirling_tuple(p+1)
    rule = :(_Jet(v, od, nc, iv, g, col, weights, marks, $p, $l, $stirling))

    return EquivariantClass( rule, eval( :((v, od, nc, iv, g, col, weights, marks) -> $rule )))
end

function _Jet(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, null::Int64, null2::Int64, null3::Int64, null4::Int64, null5::Int64, p::Int64, q::ToricLineBundle, st::Tuple{Vararg{Int64}})::Cycle
    
    try 
        1 > n_marks && error("Jet is defined when there is at least one mark")
        0 > p       && error(string("exponents of Jet must be nonnegative, correct ", p))
    catch e
        printstyled(stderr,"ERROR: ", bold=true, color=:red)
        printstyled(stderr,sprint(showerror,e), color=:light_red)
        println(stderr)
        return error_cycle()
    end
    psi_deg = p==0 ? 0 : 1
    return Cycle(p+1, psi_deg)
end

function _Jet(v::NormalToricVariety, od::Dict{Tuple{Int64, Int64}, T}, nc::Dict{Int64, Vector{Int64}}, iv::Dict{Tuple{Int64,Int64},CohomologyClass}, g::Graph{Undirected}, col::Tuple{Vararg{Int64}}, weights::Tuple{Vararg{Int64}}, marks::Tuple{Vararg{Int64}}, p::Int64, l::ToricLineBundle, st::Tuple{Vararg{Int64}})::T
    
    local s1 = zero(od[first(keys(od))])

    e = _ev(v, od, nc, iv, g, col, weights, marks, 1, l)

    for h in 0:p
        s1 += st[p+1-h]*(e^(p+1-h))*_Psi(v, od, nc, iv, g, col, weights, marks, [h])
    end

    return s1
end