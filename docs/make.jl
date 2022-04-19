using Documenter
using DipoleLattice

DocMeta.setdocmeta!(DipoleLattice, :DocTestSetup, :(using DipoleLattice); recursive=true)

makedocs(
    modules=[DipoleLattice],
    authors="Vahagn Mkhitaryan",
    repo="https://github.com/vmkhit/DipoleLattice.jl.git",
    sitename="DipoleLattice.jl",
    source = "src",
    build = "build",
    clean = true,
    doctest = true,
    highlightsig=true,
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Dipole Lattice" => "index.md",
        "Reflection and Transmission"=>"rt.md",
        "Lattice"=>"Lattice.md",
        "Multicomponent"=> "multilattice.md",
        "Supported Lattice"=> "supportedlattice.md",
        "TMM" => "tmm.md",
        "Lattice Dispersion" => "dispersion.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/vmkhit/DipoleLattice.jl.git"
)
