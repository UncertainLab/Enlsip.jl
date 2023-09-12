push!(LOAD_PATH,"../src/")

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
    versions = ["stable" => "v^", "v0.9.0"],
    devbranch = "master",
    tag_prefix = "v0.9.0"
)
