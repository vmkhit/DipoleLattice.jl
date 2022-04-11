
include("SphericalBessel.jl")

function Dn(n, z)
  return RiccatiPsip(n, z)/RiccatiPsi(n, z)
end

function sp_ab(l , λ, ɛμh, ɛμs, rs)
  m = sqrt(prod(ɛμs))/sqrt(prod(ɛμh))
  k = 2.0*pi*sqrt(prod(ɛμh))/λ
  x = k*rs
  al = (ɛμh[2]*(m^2)*sp_besselj(l, m*x)*RiccatiPsip(l, x) - ɛμs[2]*sp_besselj(l, x)*RiccatiPsip(l,m*x))/
       (ɛμh[2]*(m^2)*sp_besselj(l, m*x)*RiccatiXip(l, x) - ɛμs[2]*sp_hankelh1(l, x)*RiccatiPsip(l,m*x))

  bl = (ɛμs[2]*sp_besselj(l, m*x)*RiccatiPsip(l, x)-ɛμh[2]*sp_besselj(l, x)*RiccatiPsip(l, m*x))/
       (ɛμs[2]*sp_besselj(l, m*x)*RiccatiXip(l, x)- ɛμh[2]*sp_hankelh1(l, x)*RiccatiPsip(l, m*x))

  return al, bl
end

function Sigma_Sphere(λ, ɛμh, ɛμs, rs, id = 0)
  @assert id == 0 || id == 1
  #id = 1 calculates only the  cross-sections
  #     of dipole mode
  kh = 2.0*pi*sqrt(prod(ɛμh))/λ
  cft = 2.0*pi/kh/kh
  x = kh*rs
  if id == 1
    Nmax = 1
  else
    Nmax = round(Int, x + 4.05*x^(1/3) + 2)
  end
  sig_ext = 0.0
  sig_scat = 0.0
  for n = 1:Nmax
    an, bn = sp_ab(n, λ, ɛμh, ɛμs, rs)
    mult = (2n+1)
    sig_scat += mult*(abs2(an) + abs2(bn))
    sig_ext  += mult*real(an+bn)
  end
  return cft*sig_scat, cft*sig_ext
end

function sp_coated_ab(l, λ, ɛ, μ, a)
  # ɛ = [ɛ_sphere, ɛ_shell, ɛ_host]
  # a = [r_sphere, r_shell]

  n1 = sqrt(μ[1]*ɛ[1])
  n2 = sqrt(μ[2]*ɛ[2])
  nh = sqrt(μ[3]*ɛ[3])
  k = 2.0*pi/λ
  x = a[1]*k
  y = a[2]*k

  Al = (n1*RiccatiPsi(l,n1*x)*RiccatiPsip(l,n2*x)-n2*RiccatiPsi(l, n2*x)*RiccatiPsip(l, n1*x))/
       (n1*RiccatiPsi(l,n1*x)*RiccatiChip(l,n2*x)-n2*RiccatiChi(l, n2*x)*RiccatiPsip(l, n1*x))

  Bl = (n2*RiccatiPsi(l,n1*x)*RiccatiPsip(l,n2*x)-n1*RiccatiPsi(l, n2*x)*RiccatiPsip(l, n1*x))/
       (n2*RiccatiPsi(l,n1*x)*RiccatiChip(l,n2*x)-n1*RiccatiChi(l, n2*x)*RiccatiPsip(l, n1*x))

  al = (nh*RiccatiPsi(l,nh*y)*(RiccatiPsip(l, n2*y)-Al*RiccatiChip(l, n2*y))-n2*RiccatiPsip(l, nh*y)*(RiccatiPsi(l, n2*y)-Al*RiccatiChi(l, n2*y)))/
       (nh*RiccatiXi(l,nh*y)*(RiccatiPsip(l, n2*y)-Al*RiccatiChip(l, n2*y))-n2*RiccatiXip(l, nh*y)*(RiccatiPsi(l, n2*y)-Al*RiccatiChi(l, n2*y)))

  bl = (n2*RiccatiPsi(l,nh*y)*(RiccatiPsip(l, n2*y)-Bl*RiccatiChip(l, n2*y))-nh*RiccatiPsip(l, nh*y)*(RiccatiPsi(l, n2*y)-Bl*RiccatiChi(l, n2*y)))/
       (n2*RiccatiXi(l,nh*y)*(RiccatiPsip(l, n2*y)-Bl*RiccatiChip(l, n2*y))-nh*RiccatiXip(l, nh*y)*(RiccatiPsi(l, n2*y)-Bl*RiccatiChi(l, n2*y)))

  return al, bl

end

function sigma_sp_coated(λ, ɛ, μ, a, id = 1)
  kh = 2.0*pi/λ
  cft = 2.0*pi/kh/kh
  x = kh*a[end]
  if id == 1
    Nmax = 1
  else
    Nmax = round(Int, x + 4.05*x^(1/3) + 2)
  end
  sig_ext = 0.0
  sig_scat = 0.0
  for n = 1:Nmax
    an, bn = sp_coated_ab(n, λ, ɛ, μ, a)
    mult = (2n+1)
    sig_scat += mult*(abs2(an) + abs2(bn))
    sig_ext  += mult*real(an+bn)
  end
  return cft*sig_scat, cft*sig_ext
end



function sp_coated_ab_new(n, λ, ɛ, μ, a)
  m1 = sqrt(μ[1]*ɛ[1])/sqrt(μ[3]*ɛ[3])
  m2 = sqrt(μ[2]*ɛ[2])/sqrt(μ[3]*ɛ[3])
  m = m2/m1
  k = 2.0*pi/λ
  x = a[1]*k
  y = a[2]*k
  An = RiccatiPsi(n, m2*x)*(m*Dn(n, m1*x)-Dn(m2*x))/(m*Dn(m1*x)*RiccatiChi(n,m2*x)- RiccatiChip(n,m2*x))
  Bn = RiccatiPsi(n, m2*x)*(Dn(n, m1*x)/m-Dn(m2*x))/(Dn(m1*x)*RiccatiChi(n,m2*x)/m- RiccatiChip(n,m2*x))
  Dnt = (Dn(n, m2*y)-An*RiccatiChip(n,m2*y)/RiccatiPsi(n, m2*y))/(1-An*RiccatiChi(n,m2*y)/RiccatiPsi(n, m2*y))
  Gnt = (Dn(n, m2*y)-Bn*RiccatiChip(n,m2*y)/RiccatiPsi(n, m2*y))/(1-Bn*RiccatiChi(n,m2*y)/RiccatiPsi(n, m2*y))
  an = ((Dnt/m2+n/y)*RiccatiPsi(n, y)-RiccatiPsi(n-1, y))/((Dnt/m2+n/y)*RiccatiXi(n, y)-RiccatiXi(n-1, y))
  bn = ((Gnt*m2+n/y)*RiccatiPsi(n, y)-RiccatiPsi(n-1, y))/((Gnt*m2+n/y)*RiccatiXi(n, y)-RiccatiXi(n-1, y))
  return an, bn
end


##################### Cylinder #################################

function cyl_ab(l, ϕ, λ, ɛμh, ɛμc, rc)
  x= 2.0*pi*rc/λ
  m = sqrt(prod(ɛμc)/prod(ɛμh))
  ξ = x*sin(ϕ)
  η = x*sqrt(m^2-cos(ξ)^2+0.0im)
  D = l*cos(ϕ)*η*besselj(l,η)*hankelh1(l, ξ)*(ξ^2/η^2-1.0)
  B = ξ*(m^2 *ξ*besseljp(l,η)*besselj(l,ξ)-η*besselj(l,η)*besseljp(l, ξ))
  C = l*cos(ϕ)*η*besselj(l, η)*besselj(l,ξ)*(ξ^2/η^2-1.0)
  V = ξ*(m^2 *ξ*besseljp(l,η)*hankelh1(l,ξ)-η*besselj(l, η)*hankelh1p(l,ξ))
  W = 1.0im*ξ*(η*besselj(l,η)*hankelh1p(l,ξ)-ξ*besseljp(l,η)*hankelh1(l, ξ))
  U = (W*V+1.0im*D^2)
  al = (C*V-B*D)/U
  bl = (W*B+1.0im*D*C)/U
  return al, bl
end

############# dipole Polarizabilities  ###############
function αe_stat(λ, ϵh, ϵp, r)
    k = 2.0*pi/λ
    α0=(ϵp-ϵh)*r^3/(ϵp + 2.0*ϵh)
    return 1.0/(1.0/α0 - 2.0im*k^3/3.0)
end

function αm_stat(λ, ϵh, ϵp, r)
    k = 2.0*pi/λ
    α0 = ((k*r)^2/(30.0))*(ϵp/ϵh-1.0)*r^3
    return 1.0/(1.0/α0 - 2.0im*k^3/3.0)
end

function alpha_sp(λ, ɛμh, ɛμs, rs, id = 0)
  @assert id == 0 || id ==1 || id ==2
  # id = 0 returns: ale, elm
  # id = 1 returns 6x6 alem
  # id = 2, returns inv(alem)
  # ɛμh = [ɛ, μ ] for host medium
  # ɛμs = [ɛ, μ ] for sphere
  kh = 2.0*pi*sqrt(prod(ɛμh))/λ
  cft = 1.5im/kh^3
  a1, b1 = sp_ab(1 ,λ, ɛμh, ɛμs, rs)
  if id == 0
    return cft*a1, cft*b1
  elseif id == 1
    return [cft*a1.*eye(3) 1e-15*eye(3); 1e-15*eye(3) cft*b1.*eye(3)]
  else
    return [eye(3)./(cft*a1) 1e-15*eye(3); 1e-15*eye(3) eye(3)./(cft*b1)]
  end
end


function alpha_sp_coated( λ, ɛ, μ, a, id = 0)
  @assert id == 0 || id ==1 || id ==2
  kh = 2.0*pi*sqrt(ɛ[3]*μ[3])/λ
  cft = 1.5im/kh^3
  a1, b1 = sp_coated_ab(1, λ, ɛ, μ, a)
  if id == 0
    return cft*a1, cft*b1
  elseif id == 1
    return [cft*a1.*eye(3) zeros(3, 3); zeros(3,3) cft*b1.*eye(3)]
  else
    return [eye(3)./(cft*a1) zeros(3, 3); zeros(3,3) eye(3)./(cft*b1)]
  end
end



#################### Cylinder #################
function alpha_cy(λ, ϕ, ɛμh, ɛμc, rc)
  kh = 2.0*pi*sqrt(prod(ɛμh))/λ
  cft = 1.5im/kh^3
  a1, b1 = cyl_ab(1, ϕ, λ, ɛμh, ɛμc, rc)
  return cft*a1, cft*b1
end


########### small Ellipsoid, electrostatic approx #################
function Depfactor(a,b,c)
    L= Array(Float64, 3)
  for (i, xi)=enumerate([a, b, c])
    L[i]=0.5*(a*b*c)*quadgk(q-> 1.0/((xi^2+q)*sqrt((a^2+q)*(b^2+q)*(c^2+q))),0.0, Inf)
  end
  return L
end

function alpha_ellips(λ, ϵh, ϵp, a,b,c)
  V = 4.0*pi*a*b*c/3.0
  β = ϵh/(ϵp-ϵh)
  L = Depfactor(a,b,c)
  al = [V/(β+L[i]) for i=1:3]
  return diagm(al)
end


#=
###############            TEST                  ####################
const rs = 100.0
const ɛh = 1.0
const μh = 1.0
const μs = 1.0
const ɛs = 25.0

lam = linspace(200.0, 2000.0, 3000)
al_e = Array(Complex{Float64}, length(lam))
al_m = Array(Complex{Float64}, length(lam))
al_e_app = Array(Complex{Float64}, length(lam))
al_m_app = Array(Complex{Float64}, length(lam))


al_ec = Array(Complex{Float64}, length(lam))
al_mc = Array(Complex{Float64}, length(lam))
al_cy_e = Array(Complex{Float64}, length(lam))
al_cy_m = Array(Complex{Float64}, length(lam))

for (i, λ) = enumerate(lam)
  al_e[i], al_m[i] = alpha_sp(λ, [ɛh, μh], [ɛs, μs], rs)

  #al_ec[i] = alpha_sp_coated(λ, [ɛs, ɛh, ɛh], ones(3), [rs, rs+20.0])[1]
  #al_mc[i] = alpha_sp_coated(λ, [ɛs, ɛh, ɛh], ones(3), [rs, rs+20.0])[2]
  #al_cy_e[i] = alpha_cy(λ, pi/2, [ɛh, μh], [ɛs, μs], rs)[1]
  #al_cy_m[i] = alpha_cy(λ, pi/2, [ɛh, μh], [ɛs, μs], rs)[2]
end

using PyPlot
plot(lam, real(al_e), label = L"\Re{(\alpha_{ee}})")
plot(lam, imag(al_e), label = L"\Im{(\alpha_{ee}})")
plot(lam, real(al_m), label = L"\Re{(\alpha_{mm}})")
plot(lam, imag(al_m), label = L"\Im{(\alpha_{mm}})")


#plot(lam, real(al_e_app), label = L"\Re{(\alpha_{ee}^{app}})")
#plot(lam, imag(al_e_app), label = L"\Im{(\alpha_{ee}^{app}})")
#plot(lam, real(al_m_app), label = L"\Re{(\alpha_{mm}^{app}})")
#plot(lam, imag(al_m_app), label = L"\Im{(\alpha_{mm}^{app}})")

#plot(lam, real(al_ec), label = L"\Re{(\alpha_{ee}^{coated}})")
#plot(lam, imag(al_ec), label = L"\Im{(\alpha_{ee}^{coated}})")
#plot(lam, real(al_mc), label = L"\Re{(\alpha_{mm}^{coated}})")
#plot(lam, imag(al_mc), label = L"\Im{(\alpha_{mm}^{coated}})")

#plot(lam, real(al_cy_e), label = L"\Re{(\alpha_{ee}^{cy}})")
#plot(lam, imag(al_cy_e), label = L"\Im{(\alpha_{ee}^{cy}})")
#plot(lam, real(al_cy_m), label = L"\Re{(\alpha_{mm}^{cy}})")
#plot(lam, imag(al_cy_m), label = L"\Im{(\alpha_{mm}^{cy}})")

legend()

=#
