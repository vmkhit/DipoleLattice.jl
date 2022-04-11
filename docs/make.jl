using Documenter
using DipoleLattice

makedocs(
    sitename = "DipoleLattice",
    authors="Vahagn Mkhitaryan",
    modules = [DipoleLattice]
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/vmkhit/DipoleLattice.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Lattice"=>"lattice.md",
        "Multicomponent"=> "multilattice.md",
        "Supported Lattice"=> "supportedlattice.md",
        "TMM" => "tmm.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/vmkhit/DipoleLattice.jl.git"
)
