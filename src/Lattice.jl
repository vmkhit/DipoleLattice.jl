"""
    lattice2D
  This data struct constructs 2D lattice object and outputs its properties.
  # Arguments
  - `id::String`: String that specify one of the five 2D bravis lattices or a gneral 2D lattice as described below.
"""
struct Lattice2D
  id::String
  R::Array{<:Real, 2}
  G::Array{<:Real, 2}
  AreaWZC::Real
  AreaBZ::Real
  function  Lattice2D(id::String, args...)
    if id == "S"    # Square
      a = args[1]
      R = [a 0; 0 a]
      AreaWZC = a*a
    elseif id == "H" # hexagonal
      a = args[1]
      R = [0 -a; sqrt(3)*a/2 a/2]
      AreaWZC = a*a*sqrt(3)/2
    elseif id == "RC"        # rectangular centered
      a, b = args[1:2]
      R = [a/2 b/2; -a/2 b/2]
      AreaWZC = a*b
    elseif id == "RP"         # Rectangular primitive
      a, b = args[1:2]
      R = [a 0; 0 b]
      AreaWZC = a*b
    elseif id == "O"    # Oblique
      a, b, α = args[1:3]
      R = [a 0; b*cos(α) b*sin(α)]
      AreaWZC = a*b*sin(α)
    elseif isempty(id)  # lattice is specified with unit vectors
      a1  = args[1]
      a2  = args[2]
      R = [a1'; a2']
      AreaWZC = abs(a1[1]*a2[2] - a1[2]*a2[1])
    else
      throw("Lattice type not defined")
    end
    g1 = 2π*[ R[2, 2], -R[2, 1]] /( R[1, 1]*R[2, 2] - R[1, 2]*R[2, 1])
    g2 = 2π*[-R[1, 2],  R[1, 1]] /(-R[2, 1]*R[1, 2] + R[2, 2]*R[1, 1])
    G = [g1'; g2']
    AreaBZ = 4π^2/AreaWZC
    return new(id, R, G, AreaWZC, AreaBZ)
  end
end

# High-Throughput Computational Screening of Two-Dimensional Semiconductors
# Symmetry-restricted phase transitions in two-dimensional solids
"""
  ConstructWZC(R::Array{<:Real})
# Arguments
- `R::Array{<:Real}`: Unit cell vectors of the 2D lattice, a1= R[1, :], a2 = R[2, :]
This function implements the construction of the BZ of a 2D lattice following the procedure provided in reference
Thompson, I., and Linton, C. M. (2010). "Guided surface waves on one-and two-dimensional arrays of spheres".,
SIAM Journal on Applied Mathematics, 70(8), 2975-2995.
"""
function ConstructWZC(R::Array{<:Real})

  a = sqrt(R[1, 1]^2 + R[1, 2]^2)
  b = sqrt(R[2, 1]^2 + R[2, 2]^2)
  s1_dot_s2 = R[1, 1]*R[2, 1] + R[1, 2]*R[2, 2]
  ψ = asin(s1_dot_s2/(a*b))

  cψ, sψ = cos(ψ), sin(abs(ψ))

  # construct the perpendicular vectors and check the directions
  sp1 = [-R[1, 2], R[1, 1]]./a
  sp2 = [-R[2, 2], R[2, 1]]./b

  if (sp2[1]*R[1, 1] + sp2[2]*R[1, 2]) < 0
    sp2 = -sp2
  end
  if (sp1[1]*R[2, 1] + sp1[2]*R[2, 2]) < 0
    sp2 = -sp1
  end

  ss1 = sp2/(a*cψ)
  ss2 = sp1/(b*cψ)

  ns1 = sqrt(ss1[1]^2 + ss1[2]^2)
  ns2 = sqrt(ss2[1]^2 + ss2[2]^2)

  g = π*(ns1 - ns2*sψ)/cψ
  h = π*(ns2 - ns1*sψ)/cψ

  bz = zeros(6, 2)
  bz[1, :] = π.*ss2 + g.*R[1, :]./a
  bz[2, :] = π.*ss2 - g.*R[1, :]./a
  bz[3, :] = π.*ss1 + sign(ψ)*h.*R[2,:]./b
  bz[4, :] = -bz[1, :]
  bz[5, :] = -bz[2, :]
  bz[6, :] = -bz[3, :]
  # Removing the identical points and sorting with respect to angle
  angles = round.(atan.(bz[:, 2], bz[:, 1]), digits = 6)
  p = sortperm(angles)
  bz = bz[p, :]
  iu = indexin(unique(angles[p]), angles[p]) # unique angle indices
  return bz[iu, :]
end

"""
  make_k_path(verts::Array{<:Real, 2}, res::Union{Integer, Vector{<:Integer}}; close::Bool=true)
Generates linear interpolation paths between the vertices provided in the list vertex.
# Arguments
- `verts::Array{<:Real, 2}`:  Is a 2D arrar of points ``v_{ix}`` = verts[i, 1] and `v_{iy}`` = verts[i, 2].
- `res::Union{Integer, Vector{<:Integer}}`: Number of sampling points between vertices. If it is specified as an integer ``N``,
it will take a ``N`` sampling point between any pair of vertives. One can specify also list of integers with number
of points between each consequitive vertices.
- `close::Bool=true`: Boolian that defines if the path is closed or not.
"""
function make_k_path(verts::Array{<:Real, 2}, res::Union{Integer, Vector{<:Integer}}; close::Bool=true)
    nv = size(verts, 1)
    p = reshape([], 0, 2)
    newvert = zeros(2)
    if typeof(res) <: Integer
        if close
            for i = 1:nv
                ip1 = mod(i, nv) + 1
                for j = 1:res
                    newvert[1] = verts[i, 1] + (j-1)*(verts[ip1, 1] - verts[i, 1])/res
                    newvert[2] = verts[i, 2] + (j-1)*(verts[ip1, 2] - verts[i, 2])/res
                    p = [p; newvert']
                end
            end
        else
            for  i= 1:(nv-1)
                for j = 1:res
                    newvert[1] = verts[i, 1] + (j-1)*(verts[i+1, 1] - verts[i, 1])/res
                    newvert[2] = verts[i, 2] + (j-1)*(verts[i+1, 2] - verts[i, 2])/res
                    p = [p; newvert']
                end
            end
            p = [p; verts[end, :]']  # the last vertex is not attended in the loop, so we add it here
        end
    else
        if close
            @assert length(res) == nv
            for i = 1:nv
                ip1 = mod(i, nv) + 1
                for j = 1:res[i]
                    newvert[1] = verts[i, 1] + (j-1)*(verts[ip1, 1] - verts[i, 1])/res[i]
                    newvert[2] = verts[i, 2] + (j-1)*(verts[ip1, 2] - verts[i, 2])/res[i]
                    p = [p; newvert']
                end
            end
        else
            @assert length(res) == (nv-1)
            for i = 1:(nv-1)
                for j = 1:res[i]
                    newvert[1] = verts[i, 1] + (j-1)*(verts[i+1, 1] - verts[i, 1])/res[i]
                    newvert[2] = verts[i, 2] + (j-1)*(verts[i+1, 2] - verts[i, 2])/res[i]
                    p = [p; newvert']
                end
            end
            p = [p; verts[end, :]']
        end
    end
    return p
end

"""
  grid_in_polygon(v::Array{<:Real, 2}, n::Integer)
Generate points in a polygon. The algorithm devide the polygon into triangles by connecting the centroid of the polygon to each vertex.
Then each triangle is subdevided into smaller triangles by taking `n` points along each edge.
# Arguments
- `v::Array{<:Real, 2}`: Vertices of the polygon, specified as ``v_{ix}`` = verts[i, 1] and `v_{iy}`` = verts[i, 2].
- `n::Integer`: Number of subdivision along each edge of sub-triangles.
"""
function grid_in_polygon(v::Array{<:Real, 2}, n::Integer)
    ## POLYGON_GRID _POINTS computes points on a polygonal grid.
    #  Parameters:
    #    Input, real v[nv,v2], the coordinates of the vertices.
    #    Input, integer n, the number of subintervals.
    nv = size(v, 1)
    # Determine the centroid, and use it as the first grid point.
    vc = zeros(1, 2)
    vc[1, 1] = sum(v[1:nv, 1])/nv
    vc[1, 2] = sum(v[1:nv, 2])/nv
    # Initiate grid container
    xyg = reshape([], 0, 2)

    xyg = [xyg; vc]

    # Consider each triangle formed by two consecutive vertices and the centroid,
    # but skip the first line of points.
    for l = 1:nv
        lp1 = mod(l, nv) + 1
        for i = 1:n
            for j = 0:(n - i)
                k = n - i - j
                newv = (i*v[l, :] + j*v[lp1, :] + k*vc[1, :])/n
                xyg = [xyg; newv']
            end
        end
    end
    return xyg
end

"""
  PlotBZ(L::Lattice2D)
# Arguments
  - `L::Lattice2D`: 2D lattice object
"""
function PlotBZ(L::Lattice2D)
  bz = ConstructWZC(L.R)
  Gv = L.G
  #g1 = sqrt(Gv[1, 1]^2 + Gv[1, 2]^2)
  #g2 = sqrt(Gv[2, 1]^2 + Gv[2, 2]^2)
  n1 = [-Gv[1, 2], Gv[1, 1]]
  n2 = [-Gv[2, 2], Gv[2, 1]]

  pos1 = 0.7.*Gv[1, :] + 0.12*n1
  pos2 = 0.7.*Gv[2, :] + 0.12*n2

  fig, ax = PyPlot.subplots(1, 1)
  ax[:fill](bz[:, 1], bz[:, 2], facecolor="none", edgecolor="purple", linewidth=3)

  ax[:set_aspect]("equal")
  ax[:arrow](0, 0, Gv[1, 1], Gv[1, 2], fc="k", ec="k", head_width=0.3, head_length=0.8, width=0.05, overhang=0.1)
  ax[:arrow](0, 0, Gv[2, 1], Gv[2, 2], fc="k", ec="k", head_width=0.3, head_length=0.8, width=0.05, overhang=0.1)

  ax[:text](pos1[1], pos1[2] , "g\$_1\$", fontsize = 16)
  ax[:text](pos2[1], pos2[2], "g\$_2\$", fontsize = 16)

  ax[:tick_params](color="none", labelcolor="none")
  #ax.spines["top"][:set_edgecolor]("none")
  #ax.spines["bottom"][:set_edgecolor]("none")
  #ax.spines["left"][:set_edgecolor]("none")
  #ax.spines["right"][:set_edgecolor]("none")
  PyPlot.tight_layout()
end
