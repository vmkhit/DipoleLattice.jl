include("tmat.jl")
include("Mie.jl")
include("materials.jl")


function CircleIntersect(c1::Vector{<:Real}, R1::Real, c2::Vector{<:Real}, R2::Real)
  # res == 0   : "Circles touch each other".
  # res == 1   : "Circles do not touch each other."
  # res == -1  : "Circles intersect"
  distSq = (c1[1] - c2[1])^2 + (c1[2] - c2[2])^2
  RR = (R1 + R2) * (R1 + R2)
  if distSq == RR
    res  = 0
  elseif distSq > RR
    res =  1
  else
    res = -1
  end
  return res
end

function SphereIntersect(c1::Vector{<:Real}, R1::Real, c2::Vector{<:Real}, R2::Real)
  # res == 0   : "Circles touch each other".
  # res == 1   : "Circles do not touch each other."
  # res == -1  : "Circles intersect"
  distSq = (c1[1] - c2[1])^2 + (c1[2] - c2[2])^2 + (c1[3] - c2[3])^2
  RR = (R1 + R2) * (R1 + R2)
  if distSq == RR
    res  = 0
  elseif distSq > RR
    res =  1
  else
    res = -1
  end
  return res
end

function CheckClusterIntersect(Centers, Radii)
  ns = length(Radii)
  for i = 1:ns-1
    for j = (i+1):ns
      t = SphereIntersect(Centers[i, :], Radii[i], Centers[j, :], Radii[j])
      if t < 0
        throw("Spheres $i and $j intersect")
        break
      end
    end
  end
end


function PlotSpheres(Centers, Radii; nuv =  100)

  ns = length(Radii)
  u = range(0, 2π, length = nuv)
  v = range(0, π, length = nuv)

  x = zeros(ns, nuv, nuv)
  y = zeros(ns, nuv, nuv)
  z = zeros(ns, nuv, nuv)

  p = 0
  for (i, ui) = enumerate(u)
    for (j, vj) = enumerate(v)
      p = p + 1
      for l = 1:ns
        x[l, i, j] = Radii[l]*cos(ui)*sin(vj) - Centers[l, 1]
        y[l, i, j] = Radii[l]*sin(ui)*sin(vj) - Centers[l, 2]
        z[l, i, j] = Radii[l]*cos(vj) - Centers[l, 3]
      end
    end
  end

  fig = figure()
  ax = fig[:add_subplot](projection = "3d")
  cl  = pyimport("matplotlib.colors")
  ls =  cl.LightSource(45, 45)
  ax.set_aspect("auto")
  ax.set_box_aspect((1, 1, 0.4))

  # Shade data, creating an rgb array
  for i = 1:ns
    rgb = ls.shade(z[i, :, :], ColorMap("Reds"))
    ax[:plot_surface](x[i, :, :], y[i, :, :], z[i, :, :], rstride=1, cstride=1, linewidth=0,
                                                          antialiased=true, facecolor="black")
  end
end


function alpha_disc(λ::Real, a::Real, h::Real,  εμm::Vector{<:Number}, εμh::Vector{<:Number})
    kh = 2*π*sqrt(εμh[1]*εμh[2])/λ
    L, R = 2*a, 2*a/h
    e1 = -0.479 - 1.36*R^(0.872)
    V1V = 0.944
    a12 = 7.05/(1 - e1)
    a14 = -10.9/R^(0.98)
    VL3 = pi*(4 + 3*(R - 1)*(2*R+pi - 2))/(24*R^3)
    V1 = V1V*VL3*L^3
    s = L*kh/(2*pi)
    As = a12*s^2 + (4.0im*pi^2)*V1*s^3/(3*L^3) + a14*s^4
    alpha = (εμm[1]*V1/4π)/(1.0/(εμm[1]/εμh[1] - 1) - 1.0/(e1 - 1) - As)
  return alpha
end


function eps_Drude(λ, ϵb, ωp, γ)
  ω =  1234.94193/λ
  ϵ = ϵb - ωp^2/(ω*(ω + im*γ))
  return ϵ
end
