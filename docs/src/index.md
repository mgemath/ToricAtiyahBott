# ToricAtiyahBott.jl


```@contents
```

```@docs
ToricAtiyahBott
```
In order to install the package, type:
```julia-repl
julia> using Pkg
julia> Pkg.add("ToricAtiyahBott")
```
After the installation, simply type:
```julia-repl
julia> using ToricAtiyahBott
```
every time you want to use the program.
This package requires Oscar, so make sure that you can use Oscar before installing this package. See https://www.oscar-system.org/install/.

## The function IntegrateAB
This is the main function of the package.
```@docs
IntegrateAB
```

## Equivariant Classes
Here we list all equivariant classes currently supported by the package.
```@docs
ev
Psi
push_ev
R1_ev
Jet
class_one
push_omega
```

## Other Functions
```@docs
vir_dim_M
codimension
is_zero_cycle
a_point
moment_graph
```

## Index

```@index
```