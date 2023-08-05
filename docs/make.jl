using Documenter
using ToricAtiyahBott

makedocs(
    sitename = "ToricAtiyahBott",
    format = Documenter.HTML(),
    modules = [ToricAtiyahBott],
    doctest = false,
)

deploydocs(
    repo = "github.com/mgemath/ToricAtiyahBott.jl.git",
    versions = nothing,
    #=,
    target = "build",
    push_preview = true,=#
)