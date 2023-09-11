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
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/UncertainLab/Enlsip.jl.git",
    devbranch = "master"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
