using Documenter
using Enlsip

makedocs(
    sitename = "Enlsip.jl",
    format = Documenter.HTML(),
    modules = [Enlsip],
    pages = [
        "Home" => "index.md",
        "Usage" => "tutorial.md",
        "Reference" => "reference.md"
    ]
)

deploydocs(
    repo = "github.com/UncertainLab/Enlsip.jl.git",
    devbranch = "master"
)
