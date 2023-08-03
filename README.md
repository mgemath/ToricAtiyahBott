# AtiyahBott.jl
[![Doc](https://img.shields.io/badge/docs-stable-blue.svg)](https://mgemath.github.io/ToricAtiyahBott.jl/)

This package is a work-in-progress implementation of the Atiyah-Bott residue formula for the moduli space of genus 0 stable maps to a toric variety in the Julia language.<br>
Full documentation is available here: https://mgemath.github.io/ToricAtiyahBott.jl/.

## Installation
This package requires Oscar, so make sure that you can use Oscar before installing this package. See https://www.oscar-system.org/install/.
In order to install this package, type:
```julia-repl
julia> using Pkg
julia> Pkg.add("ToricAtiyahBott")
```
After the installation, simply type:
```julia-repl
julia> using ToricAtiyahBott
```
every time you want to use the program.

To use our code, you need to define the following: X (a toric variety), beta (the cohomology class of a curve), m (the number of marks), P (an equivariant class). See documentation for some examples.

The full list of the currently supported equivariant classes is the following:
```julia
ev(j, cc), ev(j, l) (pull back of the class cc (or line bundle l) with respect to the ev_j)
push_ev(l)  (push forward with respect to the forgetful map of the pull back of l)
R1_ev(l)           (first derived functor of direct image of the pull back of l)
Psi(a)        (cycle of psi-classes)
Jet(p,l)      (Euler class of the jet bundle J^p with respect to l)
```
Brief descriptions on these functions can be obtained through the standard help functionality of Julia by typing "?" and then the name of the function.
```julia-repl
help?> Psi
```
# Citing

We encourage you to cite our work if you have used our package. See "Cite this repository" on this page.