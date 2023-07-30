using Documenter
using ToricAtiyahBott

makedocs(
    sitename = "ToricAtiyahBott",
    format = Documenter.HTML(),
    modules = [ToricAtiyahBott]
)

deploydocs(
    repo = "github.com/mgemath/ToricAtiyahBott.jl.git",
    versions = "v#.#",
    #=,
    target = "build",
    push_preview = true,=#
)