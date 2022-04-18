
```@docs
  Lattice2D
```
```@raw html
<img src="../Pictures/lattice2D.png" width="75%"/>
```
# Description
  The id argument takes the following values, id="S" for square, id="H" for hexagonal,
  id="RC" for centered rectangular, id="RP" for primitive recangular and id="O" for
  oblique lattices. If id is specified as an empty string id="" one has to specify the
  direct lattice vectors ``a_{1}`` and ``a_{2}`` as an optional arguments. If id = "S"
  or "H" one has to specify the lattice vector length as an optional argument. Similarly,
  for the case when id="RC" or id="RP" one has to dpecify the lengths of the lattice vectors
  along each orthogonal direction. For Oblique lattice, besides the lattice vector lengths
  one has to specify also the angle between these vectors as the 3rd optional argument.

| Column One | Column Two | Column Three |
|:---------- | ---------- |:------------:|
| Row `1`    | Column `2` |              |
| *Row* 2    | **Row** 2  | Column ``3`` |

# Examples
* To define a **square lattice** with unit vector 5 we can call.
```@julia-repl
  julia> L = Lattice("S", 5)
```
* **Hexagonal lattice** with lattice vector length 4 can be defined as
```@julia-repl
  julia> L = Lattice("S", 4)
```
* **Rectangular lattices** with unit vector lengths `a = 2`, `b = 4` can be defined as.
```@julia-repl
  julia> L = Lattice("RC", a, b)
```
## or
```@julia-repl
  julia> L = Lattice("RP", a, b)
```
* The **oblique lattice** takes addtional argument ``\psi`` the angle between the vectors,
  assuming ``a_{1}`` is along the `x`-axis.

<span style="color:green;font-weight:700;font-size:20px">
    Hi Vahagn
</span>

```@julia-repl
  julia> L = Lattice("O", a, b, ψ)
```
* Finally the **general lattice** with unit vectors ``a_{1} = [a_{1x}, a_{1y}]`` and
``a_{1} = [a_{2x}, a_{2y}]``.
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
  grid_in_polygon(v::Array{<:Real, 2}, n::Integer)
```

```@docs
  PlotBZ(L::Lattice2D)
```
