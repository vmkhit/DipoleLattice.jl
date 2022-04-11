using Documenter
using DipoleLattice

makedocs(
    sitename = "DipoleLattice",
    authors="Vahagn Mkhitaryan",
    modules = [DipoleLattice]
    format=Documenter.HTML(),
    pages=[
        "Dipole Lattice" => "index.md",
        "Reflection and Transmission"=>"rt.md",
        "Lattice"=>"Lattice.md",
        "Multicomponent"=> "multilattice.md",
        "Supported Lattice"=> "supportedlattice.md",
        "TMM" => "tmm.md"
        "Lattice Dispersion" => "dispersion.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/vmkhit/DipoleLattice.jl.git"
)
