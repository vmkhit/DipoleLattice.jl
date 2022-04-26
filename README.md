
<div align="center">
  <img src="docs/src/Pictures/mapping.png" alt="Dipole Lattice Logo"></img>
</div>

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://vmkhit.github.io/DipoleLattice.jl/dev)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://vmkhit.github.io/DipoleLattice.jl/dev)
[![Runtests](https://github.com/vmkhit/DipoleLattice.jl/actions/workflows/Runtests.yml/badge.svg)](https://github.com/vmkhit/DipoleLattice.jl/actions/workflows/Runtests.yml)
# DipoleLattice.jl
This package implements various routines to calculate many scattering properties of 2D periodic arrays of dipole scatterers. The calculations are based on the assumption that the scttering of these particles can be described using their electric and magnetic dipole polarizabilities that are assumed to be given as an input. The package can calculate far-field reflection and transmission spectra of the array. Dispersion spectrum and the waveguide modes calculations are also implemented. The calculations can be done for suspended and supported array as well as for layers of dipole arrays and homogenous film. The codes can also calculate the decay rate and the lifetime of a single emitter located near such 2D periodic arrays of scatterers using lattice green function approach. 
