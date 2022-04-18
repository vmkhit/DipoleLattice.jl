
```@docs
  Lattice2D
```
```@raw html
<img src="../Pictures/lattice2D.png" width="75%"/>
```
# Description
  The `id` argument takes the following values. The table below shows the `id` values and corresponding lattice types.

| Lattice type              | id   |   
|:-------------------------:|:----:|
| **Square**                | "S"  |
| **Hexagonal**             | "H"  |
| **Centered Rectangular**  | "RC" |
| **Primitive Rectangular** | "RP" |
| **Oblique**               | "O"  |
| **General**               | ""   |

If `id` = "S" or "H" one has to specify the lattice vector length as an optional argument.
Similarly, for the case when `id`="RC" or `id`="RP" one has to specify the lengths of the lattice vectors
along each orthogonal direction. For Oblique lattice, besides the lattice vector lengths
one has to specify also the angle between these vectors as the 3rd optional argument.
If `id` is specified as an empty string `id`="" one has to specify the direct lattice vectors ``a_{1}`` and ``a_{2}``
as an optional arguments.

```math
f(a) = \frac{1}{2\pi}\int_{0}^{2\pi} (\alpha+R\cos(\theta))d\theta
```

# Examples
* To define a **square lattice** with unit vector length = `a = 5` we can call.
```@julia-repl
  julia> L = Lattice("S", 5)
```
* **Hexagonal lattice** with unit vector length  `a = 4` can be defined as
```@julia-repl
  julia> L = Lattice("S", 4)
```
* **Rectangular lattices** with unit vector lengths `a = 2`, `b = 4` can be defined as.
```@julia-repl
  julia> L = Lattice("RC", a, b)
```
The primitive rectangular lattice can be constructed using:
```@julia-repl
  julia> L = Lattice("RP", a, b)
```
* The **oblique lattice** takes additional argument ``\psi`` the angle between the vectors, assuming ``a_{1}`` is along the `x`-axis.

```@julia-repl
  julia> L = Lattice("O", a, b, ψ)
```
* Finally the **general lattice** with unit vectors ``a_{1} = [a_{1x}, a_{1y}]`` and ``a_{1} = [a_{2x}, a_{2y}]``.
```@julia-repl
  julia> L = Lattice("", a1, a2)
```

```@docs
  ConstructWZC(R::Array{<:Real})
```
```@docs
  make_k_path(verts::Array{<:Real, 2}, res::Union{Integer, Vector{<:Integer}}; close::Bool=true)
```

```@docs
  PlotBZ(L::Lattice2D)
```
# Example
```@julia-repl
  julia> L = Lattice("O", 1, 1, π/4), PlotBZ(L)
```
```@raw html
<img src="../Pictures/PlotBZ.png" width="75%"/>
```
