include("../multicomponent.jl")
include("../Lattice_Reflection.jl")
include("../utils.jl")

using PyPlot
pygui(true)


#==== ============ Test calculation goes here ==================================#
θ = 0
ϕ = 0
a = 1000.0

Centers = a.*[0.0 0.0 0.0; 0.7 0.7 0.0; 0.2 0.7 0.0]
Rs = a.*[0.1, 0.3, 0.2]

ϵμh = [1.0, 1.0]
mlat = lattice2D(a, a, pi/2)

# check if the spheres are crossing
CheckClusterIntersect(Centers, Rs)
PlotSpheres(Centers, Rs)

nlam = 20
λ = range(300.0, 2000.0, length = nlam)

nxy = 50
xx = range(0, a, length = nxy)

Escat_s = zeros(ComplexF64, nlam, nxy, nxy, 2, 3)
Hscat_s = zeros(ComplexF64, nlam, nxy, nxy, 2, 3)
Escat_p = zeros(ComplexF64, nlam, nxy, nxy, 2, 3)
Hscat_p = zeros(ComplexF64, nlam, nxy, nxy, 2, 3)
Einc_s = zeros(ComplexF64, nlam, nxy, nxy, 3)
Hinc_s = zeros(ComplexF64, nlam, nxy, nxy, 3)
Einc_p = zeros(ComplexF64, nlam, nxy, nxy, 3)
Hinc_p = zeros(ComplexF64, nlam, nxy, nxy, 3)

for l = 1:nlam
  k0 = 2*π/λ[l]
  kp = k0.*[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ)]
  #ϵm1 = eps_Drude(λ[l], 4.0, 9.01, 0.021)  # Ag
  #ϵm2 = eps_Drude(λ[l], 8.0, 9.17, 0.071)  # Au
  ϵm1 = 9.0
  ϵm2 = 16.0
  EpsMus = [ϵm1 1.0; ϵm2 1.0; ϵm2 1.0]
  Sp = Spheres(Centers, Rs, EpsMus)

  al_e1, al_m1 = al_em(λ[l], [Sp.EpsMu[1, 1], ϵμh[1]], [Sp.EpsMu[1, 2], ϵμh[2]], [Sp.Radii[1]])
  al_e2, al_m2 = al_em(λ[l], [Sp.EpsMu[2, 1], ϵμh[1]], [Sp.EpsMu[2, 2], ϵμh[2]], [Sp.Radii[2]])
  al_e3, al_m3 = al_em(λ[l], [Sp.EpsMu[3, 1], ϵμh[1]], [Sp.EpsMu[3, 2], ϵμh[2]], [Sp.Radii[3]])

  α1 = Make_alpha6x6(al_e1, al_m1)
  α2 = Make_alpha6x6(al_e2, al_m2)
  α3 = Make_alpha6x6(al_e3, al_m3)

  α1_ee = Make_alpha3x3(al_m1)
  α2_ee = Make_alpha3x3(al_m2)
  α3_ee = Make_alpha3x3(al_m3)

  Dp = Dipoles(Centers, [α1, α2, α3])
  for (i, xi) = enumerate(xx)
    for (j, yj) = enumerate(xx)
      Einc_s[l, i, j, :], Hinc_s[l, i, j, :], Escat_s[l, i, j, :, :], Hscat_s[l, i, j, :, :] = EHFAR_PM_DipoleCluster("s", k0, θ, ϕ, ϵμh, Dp, mlat, xi, yj, 10*a)
      Einc_p[l, i, j, :], Hinc_p[l, i, j, :], Escat_p[l, i, j, :, :], Hscat_p[l, i, j, :, :] = EHFAR_PM_DipoleCluster("p", k0, θ, ϕ, ϵμh, Dp, mlat, xi, yj, 10*a)
    end
  end
end


pcolormesh(xx, xx, real.(Escat_s[6, :, :, 1, 2]))
