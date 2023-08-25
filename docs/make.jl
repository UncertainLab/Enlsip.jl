push!(LOAD_PATH,"../src/")

using Documenter
using Enlsip

makedocs(
    sitename = "Enlsip.jl",
    format = Documenter.HTML(),
    modules = [Enlsip],
    pages = [
        "Home" => "index.md",
        "How to use" => "tutorial.md",
        "API" => "api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
