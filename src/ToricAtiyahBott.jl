"""
**ToricAtiyahBott** is a module containing an implementation of the Kontsevich-Atiyah-Bott residue formula for toric varieties in the Julia language.
The theory and the algorithm behind the package is described in the paper 
"Computations of Gromov-Witten invariants of toric varieties" by Giosuè Muratore.

"""
module ToricAtiyahBott

using Oscar
using ProgressMeter
using Combinatorics

include("Arithmetic.jl")
include("Main.jl")

include("Checks.jl")
include("Colors.jl")
include("Euler.jl")
include("ev_class.jl")
include("Jet_class.jl")
include("class_one_class.jl")

include("Marks.jl")
include("MaxCones.jl")
include("Misc.jl")
include("Multi.jl")
include("Omega.jl")
include("Psi_class.jl")
include("push_ev_class.jl")
include("push_omega_class.jl")
include("R1_ev_class.jl")
include("Rules.jl")

include("Trees.jl")

end # module ToricAtiyahBott