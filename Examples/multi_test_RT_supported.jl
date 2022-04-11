include("../src/DipoleArray_supported.jl")
include("../src/Lattice_Reflection.jl")
include("../src/utils.jl")

using PyPlot
pygui(true)


#==== ============ Test calculation goes here ==================================#
θ = pi/4
ϕ = 0
a = 1000.0

Centers = a.*[0.0 0.0 0.1; 0.75 0.3 0.1; 0.2 0.7 0.1]
Radii = a.*[0.1, 0.1, 0.1]

ϵμh = [1.0, 1.0]
mlat = lattice2D(a, a, pi/2)

# check if the spheres are crossing
CheckClusterIntersect(Centers, Radii)
PlotSpheres(Centers, Radii)

mat0 = Material()
mat1 = Material("epsmu",  3.0, 1.0)
mat2 = Material("epsmu",  2.1, 1.0)
L1 = Layer(mat0, 0)
L3 = Layer(mat1, 600)
L5 = Layer(mat1, 0)


nlam = 500
λ = range(300.0, 2000.0, length = nlam)

R_s = zeros(nlam)
T_s = zeros(nlam)
A_s = zeros(nlam)

R_p = zeros(nlam)
T_p = zeros(nlam)
A_p = zeros(nlam)


for i=1:nlam
  k0 =  2*π/λ[i]
  #ϵm1 = eps_Drude(λ[i], 4.0, 9.01, 0.021)  # Ag
  #ϵm2 = eps_Drude(λ[i], 8.0, 9.17, 0.071)  # Au
  ϵm1 = 9.0+0.2im
  ϵm2 = 4.0+0.2im
  EpsMu = [ϵm1 1.0; ϵm2 1.0; ϵm2 1.0]

  ϵm = eps_Drude(λ[i], 4.0, 9.01, 0.1)
  mat3 = Material("epsmu", ϵm, 1.0)
  L2 = Layer(mat3, 20)
  L4 = Layer(mat3, 100)
  S = Stack([L1, L2, L3, L4, L5], zeros(4))


  al_e1, al_m1 = al_em(λ[i], [EpsMu[1, 1], ϵμh[1]], [EpsMu[1, 2], ϵμh[2]], [Radii[1]])
  al_e2, al_m2 = al_em(λ[i], [EpsMu[2, 1], ϵμh[1]], [EpsMu[2, 2], ϵμh[2]], [Radii[2]])
  al_e3, al_m3 = al_em(λ[i], [EpsMu[3, 1], ϵμh[1]], [EpsMu[3, 2], ϵμh[2]], [Radii[3]])

  α1 = Make_alpha6x6(al_e1, al_m1)
  α2 = Make_alpha6x6(al_e2, al_m2)
  α3 = Make_alpha6x6(al_e3, al_m3)

  Dp = Dipoles(Centers, [α1, α2, α3])

  R_s[i], T_s[i], A_s[i] =   RT_PM_DipoleCluster_supported("s", k0, θ, ϕ, ϵμh, Dp, mlat, S)
  R_p[i], T_p[i], A_p[i] =   RT_PM_DipoleCluster_supported("p", k0, θ, ϕ, ϵμh, Dp, mlat, S)

end

begin
  fig, ax = subplots(1, 3, figsize = (15, 5))
  ax[1][:plot](λ, R_s, label = L"R_{pp}", color = "red" ,linestyle = "-")
  ax[1][:plot](λ, R_p, label = L"R_{pp}", color = "green" ,linestyle = "-")


  ax[2][:plot](λ, T_s, label = L"T_{ss}", color = "red" ,linestyle = "-")
  ax[2][:plot](λ, T_p, label = L"T_{pp}", color = "green" ,linestyle = "-")


  ax[3][:plot](λ, A_s, label = L"R_{pp}", color = "red" ,linestyle = "-")
  ax[3][:plot](λ, A_p, label = L"R_{pp}", color = "green" ,linestyle = "-")
  for i = 1:3
    #ax[i][:set_xlim]([400, 1700])
    #ax[i][:set_ylim]([-1, 2])
  end
end
