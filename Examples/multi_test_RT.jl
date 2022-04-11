using PyPlot
pygui(true)


#==== ============ Test calculation goes here ==================================#
θ = 0
ϕ = pi/4
a = 1000.0
Centers = a.*[0.0 0.0 0.0; 0.7 0.3 0.0; 0.2 0.7 0.0]
Radii = a.*[0.2, 0.2, 0.1]
ϵμh = [1.0, 1.0]  # [ε, μ]
mlat = lattice2D(a, a, pi/3)
# check if the spheres are crossing
CheckClusterIntersect(Centers, Radii)
PlotSpheres(Centers, Radii)
nlam = 500
λ = range(500.0, 3000.0, length = nlam)

R_s = zeros(nlam)
T_s = zeros(nlam)
A_s = zeros(nlam)

R_p = zeros(nlam)
T_p = zeros(nlam)
A_p = zeros(nlam)

rss = zeros(ComplexF64, nlam)
rsp = zeros(ComplexF64, nlam)
tss = zeros(ComplexF64, nlam)
tsp = zeros(ComplexF64, nlam)

rpp = zeros(ComplexF64, nlam)
rps = zeros(ComplexF64, nlam)
tpp = zeros(ComplexF64, nlam)
tps = zeros(ComplexF64, nlam)



for i = 1:nlam
  println(λ[i])

  k0 =  2*π/λ[i]

  #ϵm1 = eps_Drude(λ[i], 4.0, 9.01, 0.021)  # Ag eps_Drude(λ, ϵb, ωp[eV], γ[eV]) = ϵb- ωp^2/(ω(ω+iγ))
  #ϵm2 = eps_Drude(λ[i], 8.0, 9.17, 0.071)  # Au

  ϵm1 = 9.0
  ϵm2 = 4.0
  ϵm3 = 16.0

  EpsMu = [ϵm1 1.0; ϵm2 1.0; ϵm2 1.0]

  al_e1, al_m1 = al_em(λ[i], [EpsMu[1, 1], ϵμh[1]], [EpsMu[1, 2], ϵμh[2]], [Radii[1]])
  al_e2, al_m2 = al_em(λ[i], [EpsMu[2, 1], ϵμh[1]], [EpsMu[2, 2], ϵμh[2]], [Radii[2]])
  al_e3, al_m3 = al_em(λ[i], [EpsMu[3, 1], ϵμh[1]], [EpsMu[3, 2], ϵμh[2]], [Radii[3]])

  α1 = Make_alpha6x6(al_e1, al_m1)
  α2 = Make_alpha6x6(al_e2, al_m2)
  α3 = Make_alpha6x6(al_e3, al_m3)

  Dp = Dipoles(Centers, [α1, α2, α3])

  R_s[i], T_s[i], A_s[i] = RT_PM_DipoleCluster_diffractive("s", k0, θ, ϕ, ϵμh, Dp, mlat)
  R_p[i], T_p[i], A_p[i] = RT_PM_DipoleCluster_diffractive("p", k0, θ, ϕ, ϵμh, Dp, mlat)

  rss[i], rsp[i], tss[i], tsp[i]  = RT_PM_DipoleCluster("s", k0, θ, ϕ, ϵμh, Dp, mlat)

  rps[i], rpp[i], tps[i], tpp[i]  = RT_PM_DipoleCluster("p", k0, θ, ϕ, ϵμh, Dp, mlat)

end


begin
  fig, ax = subplots(1, 3, figsize = (15, 5))
  ax[1][:plot](λ, R_s, label = L"R_{pp}", color = "red" ,linestyle = "-")
  ax[1][:plot](λ, R_p, label = L"R_{pp}", color = "green" ,linestyle = "-")
  ax[1][:set_ylabel]("Reflectance")


  ax[2][:plot](λ, T_s, label = L"T_{ss}", color = "red" ,linestyle = "-")
  ax[2][:plot](λ, T_p, label = L"T_{pp}", color = "green" ,linestyle = "-")
  ax[2][:set_ylabel]("Transmittance")


  ax[3][:plot](λ, A_s, label = L"R_{pp}", color = "red" ,linestyle = "-")
  ax[3][:plot](λ, A_p, label = L"R_{pp}", color = "green" ,linestyle = "-")
  ax[3][:set_ylabel]("Absorbance")
  for i = 1:3
    ax[i][:set_xlabel]("Wavelength (nm)")
  end
  for i = 1:3
    #ax[i][:set_xlim]([400, 1700])
    #ax[i][:set_ylim]([-1, 2])
  end
  tight_layout()
end


begin
  fig, ax = subplots(1, 3, figsize = (15, 5))
  ax[1][:plot](λ, abs2.(rss), label = L"R_{pp}", color = "red" ,linestyle = "-")
  ax[1][:plot](λ, abs2.(rpp), label = L"R_{pp}", color = "green" ,linestyle = "-")


  ax[2][:plot](λ, abs2.(tss), label = L"T_{ss}", color = "red" ,linestyle = "-")
  ax[2][:plot](λ, abs2.(tpp), label = L"T_{pp}", color = "green" ,linestyle = "-")


  ax[3][:plot](λ, 1 .- abs2.(rss) .- abs2.(tss), label = L"R_{pp}", color = "red" ,linestyle = "-")
  ax[3][:plot](λ, 1 .- abs2.(rpp) .- abs2.(tpp), label = L"R_{pp}", color = "green" ,linestyle = "-")
  for i = 1:3
    #ax[i][:set_xlim]([400, 1700])
    #ax[i][:set_ylim]([-1, 2])
  end
end
