using LinearAlgebra
include("SphericalBessel.jl")


function hankel_plus(l, z)
  return im*sp_hankelh1(l, z)
end

function hankel_plus_p(l,z)
  return im*sp_hankelh1p(l, z)
end

#sphericalbesselj(nu,z)

function tEM_sp(l, λ, ɛ, μ , a)

  # ɛ = Vector(N+1), ɛ[N+1] = ɛh
  # μ  = Vector(N+1), μ[N+1] = μh
  # a = Vector(N)
  # array orders follow from inner sphere to outer shell

  k = 2.0*pi/λ
  ME_tot = Matrix{Float64}(I, 2, 2)
  MH_tot = Matrix{Float64}(I, 2, 2)
  for j = 1:length(a)

    ρ = k*sqrt(ɛ[j]*μ[j])*a[j]
    ρt = k*sqrt(ɛ[j+1]*μ[j+1])*a[j]

    I = sp_besselj(l,ρ)
    It = sp_besselj(l, ρt)

    H = hankel_plus(l, ρ)
    Ht = hankel_plus(l, ρt)

    U = sp_besselj(l, ρ) + ρ*sp_besseljp(l,ρ)
    V = hankel_plus(l, ρ) + ρ*hankel_plus_p(l, ρ)

    Ut = sp_besselj(l,ρt)+ρt*sp_besseljp(l,ρt)
    Vt = hankel_plus(l,ρt)+ρt*hankel_plus_p(l,ρt)

    cft_E = 1.0/(ɛ[j+1]*(It*Vt-Ut*Ht))
    cft_M = 1.0/(μ[j]*(Vt*It-Ut*Ht))

    ME_j = cft_E.*[(ɛ[j]*I*Vt - ɛ[j+1]*U*Ht) (ɛ[j]*H*Vt - ɛ[j+1]*V*Ht);
                  (-ɛ[j]*I*Ut + ɛ[j+1]*U*It) (-ɛ[j]*H*Ut + ɛ[j+1]*V*It)]

    MH_j = cft_M.*[(μ[j]*I*Vt - μ[j+1]*U*Ht) (μ[j]*H*Vt - μ[j+1]*V*Ht);
                 -(μ[j]*I*Ut - μ[j+1]*U*It) (-μ[j]*H*Ut + μ[j+1]*V*It)]

    ME_tot = ME_j*ME_tot
    MH_tot = MH_j*MH_tot
  end

  #A = 1.0/M_tot[1,1]             # this gives the field inside innermost sphere as E = A*j_{L}(R₀)
  tE = ME_tot[2,1]/ME_tot[1,1]    # this gives the field in N+1 medium as j_{L}(R_(N+1))+t*h_{L}(R_(N+1))
  tH = MH_tot[2,1]/MH_tot[1,1]

  return [tE, tH]
end

function al_em(λ, ɛ,  μ , a)
  k = 2.0*pi/λ
  return 1.5*tEM_sp(1, λ, ɛ, μ , a)/k^3
end


##################### TEST ########################

#=

lam = range(100, 3000, length = 5000)
al = Array{ComplexF64}(undef, length(lam), 2)
for (i,λ) = enumerate(lam)
  al[i, :] = al_em(λ, [12.0, 3.0, 1.0], ones(3) , [500.0, 100])
end

using PyPlot
pygui(true)

begin
  plot(lam, real(al[:, 1]), label = L"\Re(\alpha_{e})")
  plot(lam, imag(al[:, 1]), label = L"\Im(\alpha_{e})")
  plot(lam, real(al[:, 2]), label = L"\Re(\alpha_{m})")
  plot(lam, imag(al[:, 2]), label = L"\Im(\alpha_{m})")
  legend()
end

=#
