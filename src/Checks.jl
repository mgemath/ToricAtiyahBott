export 
    codimension,
    vir_dim_M,
    is_zero_cycle

import Base: *, //, /, ^, +, -, inv, one, zero

"""
    vir_dim_M(v, beta, m)

The virtual dimension of the moduli space of stable rational map to the toric variety `v`, of class `beta` with `m` marks.
# Arguments
- `v::NormalToricVariety`: the toric variety.
- `beta::CohomologyClass`: the class of the stable maps.
- `m::Int64`: the number of marks.

# Example
```julia-repl
julia> v = projective_space(NormalToricVariety, 2);
julia> beta = cohomology_class(toric_divisor(v, [1,0,0]));
julia> vir_dim_M(v,beta,0)
2
```
"""
function vir_dim_M(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64 = 0)::Int64
    
    return dim(v) - 3 - Int64(integrate(cohomology_class(canonical_divisor(v))*beta)) + n_marks
end

"""
    codimension(v, beta, m, P)

The codimension of the equivariant class `P`.
# Arguments
- `v::NormalToricVariety`: the toric variety.
- `beta::CohomologyClass`: the class of the stable maps.
- `m::Int64`: the number of marks. If omitted, it is assumed as `m=0`.
- `P`: the equivariant class.
"""
function codimension(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, P)::Int64
    
    return Base.invokelatest( P.func, v, beta, n_marks, 0, 0, 0, 0, 0).c
end

"""
    is_zero_cycle(v, beta, m, P)

Return `true` if the equivariant class `P` is a 0-cycle in the moduli space, `false` otherwise.
# Arguments
- `v::NormalToricVariety`: the toric variety.
- `beta::CohomologyClass`: the class of the stable maps.
- `m::Int64`: the number of marks.
- `P`: the equivariant class.
"""
function is_zero_cycle(v::NormalToricVariety, beta::CohomologyClass, n_marks::Int64, P_input)::Bool

    local n_results::Int64 = 0

    if isa(P_input, Array)  #we want that P is an array, possibly with only one element
        n_results = length(P_input)
    else
        n_results = 1
        #P = [P]
    end
    local P::Vector{Function} = Vector(undef, n_results)
    
    if isa(P_input, Array)  #we want that P is an array, possibly with only one element
        if typeof(P_input[1]) == EquivariantClass
            for res in eachindex(P)
                P[res] = P_input[res].func
            end
        else
            P = P_input
        end
    else
        if typeof(P_input) == EquivariantClass
            P[1] = P_input.func
        else
            P[1] = P_input
        end
    end

    for res in eachindex(P)
        local P_cycle = Base.invokelatest( P[res], v, beta, n_marks, 0, 0, 0, 0, 0)

        if !P_cycle.valid
            printstyled("Warning: ", bold=true, color=:light_yellow)
            println("the class is not valid")
            return false
        end

        if P_cycle.is_psi > 1
            printstyled("Warning: ", bold=true, color=:light_yellow)
            println("more instances of Psi has been found. Type:")
            printstyled("julia> ", bold=true, color=:light_green)
            println("?Psi")
            println("for support.")
            return false
        end

        if P_cycle.c != vir_dim_M(v, beta, n_marks)
            printstyled("Warning: ", bold=true, color=:light_yellow)
            length(P)==1 ? println("the class is not a zero cycle") : println("some classes are not zero cycles")
            return false
        end
    end
    return true
end

##############Others################

function split_cc(c::CohomologyClass)::Vector{CohomologyClass}
    p = polynomial(c)
    v = toric_variety(c)
    # ring = cohomology_ring(v)
    # if iszero(p)
    #     return [cohomology_class(v, zero(ring))]
    # end
    coeffs = [k for k in AbstractAlgebra.coefficients(p.f)]
    expos = matrix(ZZ, [k for k in AbstractAlgebra.exponent_vectors(p.f)])
    indets = gens(cohomology_ring(toric_variety(c)))
    monoms = [prod(cohomology_class(v, indets[j])^expos[k, j] for j in 1:ncols(expos)) for k in 1:nrows(expos)]
    return [coeffs[k]*monoms[k] for k in eachindex(monoms)]
end

struct Cycle #We define the structure Cycle. It keeps track of the codimension and of the fact that it contains psi-classes
    c::Int64 #The codimension. If negative, the class is not well defined
    is_psi::Int64 #this is a counter that measure how many times the class Psi is called in P. If this number is strictly greater than 1, then the program stops.
    valid::Bool #check if it is a valid cycle
end

function error_cycle()::Cycle
    return Cycle(0, 0, false)
end

function Cycle(c::Int64, is_psi::Int64)::Cycle
    return Cycle(c, is_psi, true)
end

#Cycles satisfy the following arithmetic. They are not invertible (except in codimension 0), and cycles of different codimension are not summable.
#So the Cycle(-1, 0) is a meaningless cycle.

*(P1::Cycle, P2::Cycle)::Cycle = Cycle(P1.c + P2.c, P1.is_psi + P2.is_psi, P1.valid && P2.valid)
^(P1::Cycle, n::Int64) ::Cycle = Cycle(n*P1.c, n*P1.is_psi, P1.valid) #if n=0, we got a zero cycle
+(P1::Cycle, P2::Cycle)::Cycle = P1.c == P2.c ? Cycle(P1.c, max(P1.is_psi, P2.is_psi), P1.valid && P2.valid) : error_cycle()
-(P1::Cycle, P2::Cycle)::Cycle = +(P1, P2)
*(P1::Cycle, ::Number) ::Cycle = P1
*(::Number, P1::Cycle) ::Cycle = P1
//(P1::Cycle, ::Number)::Cycle = P1
//(::Number, P1::Cycle)::Cycle = Cycle(-P1.c, -P1.is_psi, P1.valid)
/(P1::Cycle, ::Number) ::Cycle = P1
/(::Number, P1::Cycle) ::Cycle = Cycle(-P1.c, -P1.is_psi, P1.valid)
+(P1::Cycle, n::Number)::Cycle = n == 0 ? P1 : error_cycle()
-(P1::Cycle, n::Number)::Cycle = +(P1, n)
+(n::Number, P1::Cycle)::Cycle = +(P1, n)
-(n::Number, P1::Cycle)::Cycle = +(P1, n)
# +(P1::Cycle)           ::Cycle = P1
# -(P1::Cycle)           ::Cycle = P1
//(P1::Cycle,P2::Cycle)::Cycle = Cycle(P1.c-P2.c, P1.is_psi-P2.is_psi, P1.valid && P2.valid)
/(P1::Cycle, P2::Cycle)::Cycle = //(P1, P2)
one(::Cycle)           ::Cycle = Cycle(0, 0, true)
inv(P1::Cycle)         ::Cycle = Cycle(-P1.c, -P1.is_psi, P1.valid)
zero(::Cycle)          ::Int64 = 0