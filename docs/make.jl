push!(LOAD_PATH,"../src/")

using Documenter
using Enlsip

makedocs(
    sitename = "Enlsip",
    format = Documenter.HTML(),
    modules = [Enlsip]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
