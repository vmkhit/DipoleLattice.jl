using Documenter
using DipoleLattice

makedocs(
    sitename = "DipoleLattice",
    format = Documenter.HTML(),
    modules = [DipoleLattice]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
