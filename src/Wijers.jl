using SpecialFunctions
using LinearAlgebra
import GSL: sf_gamma_inc


mutable struct lattice2D
    a::Real
    b::Real
    alpha::Real
    R::Array{Real, 2}
    G::Array{Real, 2}
    AWZC::Real
    function lattice2D(a, b, alpha)
        R = zeros(2,2)
        R[1,1], R[1,2] = a, zero(a)
        R[2,1], R[2,2] = b*cos(alpha), b*sin(alpha)
        G = 2.0*pi*inv(R')
        AWZC = a*b*sin(alpha)
        return new(a, b, alpha, R, G, AWZC)
    end
end
Lattice2D(a::Real, b::Real, alpha::Real) = Lattice2D(promote(a,b,alpha)...)

const eps_lut = cat(dims = 3, [0 0 0; 0 0 1; 0 -1 0], [0 0 -1; 0 0 0; 1 0 0], [0 1 0; -1 0 0; 0 0 0]);

function ixa(a::Vector{T}, ν::Int, μ::Int) where T<:Union{Float64, Complex}
    s = zero(eltype(a))
    for l = 1:3
        s += eps_lut[μ, ν, l]*a[l]
    end
    return s
end
const dfn = [1 0 0; 0 1 0; 0 0 1];
const sfn = [0 0 0; 0 0 0; 0 0 1];
const tfn = [1 1 -1; 1 1 -1; -1 -1 1];

######################## Spatial domain Green Functions ########################


function _g0(k::Number, r::Real)
    return cis(k*r)/(4.0*pi*r)
end
function _P(z::Number)
    return (1. - 1.0/z + 1.0/z^2)
end
function _Q(z::Number)
    return (-1.0 + 3.0/z - 3.0/z^2)
end

function GE_reg(k::Number, r1::Vector{T}, r2::Vector{T}; id ="F", a::Real = 1e-3) where T<:Real
    @assert id == "F" || id == "T" || id == "L" # for full, transcerse and longitudinal parts
    ΛL = 1.5^(1/3)/a
    ΛL3 = 1.5/a^3
    ΛT = 1.5/a
    G = Array{ComplexF64}(undef, 3,3)
    r1 = r1-r2
    r = norm(r1)
    k2 = k*k

    if r < 1e-9
        for j = 1:3
            for i = 1:3
                if id == "L"
                    G[i, j] = -dfn[i, j]*ΛL3/(6.0*pi*k2)
                elseif id == "T"
                    G[i, j] = dfn[i, j]*(ΛT+im*k)/(6.0*pi)
                else
                    G[i, j] = dfn[i, j]*((ΛT+im*k)/(6.0*pi) - ΛL3/(6.0*pi*k2))
                end
            end
        end
    else
        for j = 1:3
            for i = 1:3
                if id == "L"
                    G[i, j] = -(dfn[i, j]-3.0*r1[i]*r1[j])*fL1 - r1[i]*r1[j]*fL2
                elseif id == "T"
                    G[i, j] = (dfn[i, j]-3.0*r1[i]*r1[j])*(_g0(k, r)*(_P(im*k*r)*dfn[i, j] + _Q(im*k*r)*r1[i]*r1[j]))
                            + (_P(-ΛT*r)*dfn[i, j] +_Q(-ΛT*r)*r1[i]*r1[j])*exp(-ΛT*r)/(4.0*pi*r)
                else
                    G[i, j] = -(dfn[i, j]-3.0*r1[i]*r1[j])*fL1 + r1[i]*r1[j]*fL2 +
                          (_P(-ΛT*r)*dfn[i, j] +_Q(-ΛT*r)*r1[i]*r1[j])*exp(-ΛT*r)/(4.0*pi*r) +
                          (dfn[i, j]-3.0*r1[i]*r1[j])*(_g0(k, r)*(_P(im*k*r)*dfn[i, j] + _Q(im*k*r)*r1[i]*r1[j]))
                end
            end
        end
    end
    return G
end

function _GEH(k::Number, r1::Vector{T}, r2::Vector{T}) where T<:Real
    r1 = r1-r2
    r = norm(r1)
    x = im*k*r
    g0 = cis(k*r)/(4.0im*pi*k*r^3)
    GE = GH = zeros{Complex64}(3, 3)
    for j = 1:3
        for i = 1:3
            GE[i, j] = (g0/(im*k))*((1.0 - x + x^2)*dfn[i, j] + (-3.0 + 3.0*x - x^2)*r1[i]*r1[j]/r^2)
            GH[i, j] =  g0*(-1.0 + x)*ixa(r1-r2, i, j)
        end
    end
    return GE, GH
end

######################## LATTICE SUMS ######################################

function Gammas(r::T, k::Union{T, Complex{T}}, η::T) where T<:AbstractFloat
    r2::T = r*r
    η2::T = η*η
    k2::Union{T, Complex{T}} = k*k
    kr::Union{T, Complex{T}} = k*r
    ηr::T = r*η
    EXP = cis(kr)
    EE = erfc(ηr + 0.5im*k/η)
    G1 = (4.0*ηr/sqrt(pi))*exp(-ηr*ηr + 0.25*k2/η2)
    G2 = (-1.0 + 1.0im*kr + kr*kr)*EXP
    G3 = (-3.0/r2 + 3.0im*k/r + k2)*EXP
    G4 = 3.0/r2 + 2.0*η2
    G5 = (1.0 - 1.0im*kr)*EXP
    return [EE, G1, G2, G3, G4, G5]
end

function SigmaDelta(z::T, kz::Union{T, Complex{T}}, η::T) where T<:AbstractFloat
    zη::T = z*η
    kzz::Union{T, Complex{T}} = kz*z
    kη::Union{T, Complex{T}} = 0.5*kz/η
    S = 4im*η*exp(-zη*zη + kη*kη)/sqrt(pi)
    Dp = cis(kzz)*erfc(-zη - im*kη)
    Dm = cis(-kzz)*erfc(zη - im*kη)
    return [S, Dp, Dm]
end

function spiral(n::Int)
    k::Int = ceil(Int, (sqrt(n)-1)/2)
    m::Int = (2*k+1)^2
    t::Int = 2k
    if n>=m-t  return  [k-(m-n); -k]        else m=m-t end
    if n>=m-t  return [-k; -k+(m-n)]      else m=m-t end
    if n>=m-t  return [-k+(m-n); k] else return [k; k-(m-n-t)] end
end

function Lsum(rv::Vector{T}, k::Union{T, Complex{T}}, kp::Vector{<:Real}, lattice::lattice2D; Nmax::Int = 3, self::Bool = false) where T<:Union{Float64, Int}
    Ra::Array{T} = lattice.R
    Ga::Array{T} = lattice.G
    AA::T = lattice.AWZC
    η::T = max(1.8*sqrt(pi/AA), 0.5*k/sqrt(6))
    η2::T = η*η
    k2::Union{T, Complex{T}} = k*k
    zz::T = abs(rv[3])

    #initiating the needed arrays
    Ge = Array{ComplexF32}(undef, 3,3)
    Gh = Array{ComplexF32}(undef, 3,3)
    gv = similar(kp)
    Nv = similar(rv)
    Kg = Array{Complex{T}}(undef, 3)
    cft::Complex = 1.0im*pi/AA;
    r::T = norm(rv);
    ############ 1: calculating the singular parts
    if r < 1e-15
        AL = AT = AI = Γ1 = Γ2 = Γ3 = Γ4 = Γ5 = EE = 0
        for μ = 1:3
            for ν = 1:3
                Gh[ν, μ] = complex(0.0, 0.0)
                if ν==μ
                    Ge[ν, μ] = -(4.0*η*exp(0.25*k2/η2)*(k2-η2)/(3*sqrt(pi))+(2.0im*k2*k/3.0)*erfc(-0.5im*k/η))
                    if self # to include the radiative damping due to self interaction
                        Ge[ν, μ] = Ge[ν, μ] + 2.0im*k2*k/3.0
                    end
                else
                    Ge[ν, μ] = complex(0.0, 0.0)
                end
            end
        end
    else
        EE, Γ1, Γ2, Γ3, Γ4, Γ5 = Gammas(r, k, η)
        AL = -0.5*(Γ1 + 2.0*(Γ2 - real(Γ2*EE)))/r^3
        AT = -0.5*(-Γ1*Γ4 - 2.0*(Γ3-real(Γ3*EE)))/r^3
        AI =  -0.5*(Γ1-2.0*(Γ5-real(Γ5*EE)))/r^3
        for ν = 1:3
            for μ = 1:3
                Ge[ν, μ] = AL*dfn[ν,μ]+ AT*rv[ν]*rv[μ] # minus sign is becaus the singular part is substracted
                Gh[ν, μ] = 1.0im*k*AI*ixa(rv, ν, μ)
            end
        end
    end
    ############ 2: Calculating the main sum
    N = 4*Nmax^2 + 1
    for n = 1:N
        inds = spiral(n)
        #Initiating reciprocal sum parameters
        for l = 1:2
            Kg[l] = kp[l] + inds[1].*Ga[1,l] + inds[2].*Ga[2,l]
            Nv[l] = rv[l] - inds[1].*Ra[1,l] - inds[2].*Ra[2,l]
        end
        Nv[3] = rv[3]
        Kg[3] = sqrt(k2-norm(Kg[1:2])^2+0.0im)
        Σ, Δp, Δm = SigmaDelta(zz, Kg[3], η)
        PhiG = cis(dot(Kg[1:2], rv[1:2]))
        #Initiating direct sum parameters

        r = norm(Nv);
        if abs2(inds[1])+abs2(inds[2]) == 0
            AL = AT = AI = Γ1 = Γ2 = Γ3 = Γ4 = Γ5 = EE = PhiR = 0
        else
            EE, Γ1, Γ2, Γ3, Γ4, Γ5 = Gammas(r, k, η)
            PhiR = 0.5*cis(dot(kp, Nv[1:2]))/r^3
            AL = - Γ1 + 2.0*real(Γ2*EE)
            AT = Γ1*Γ4 - 2.0*real(Γ3*EE)
            AI = - Γ1 - 2.0*real(Γ5*EE)
        end
        for μ = 1:3
            for ν = 1:3
                Ge[ν, μ] += cft*PhiG*((k2*dfn[ν, μ] - Kg[ν]*Kg[μ])*(Δp + tfn[ν, μ]*Δm)/Kg[3]  + sfn[ν,μ]*Σ) + PhiR*(AL*dfn[ν, μ] + AT*Nv[ν]*Nv[μ])
                Gh[ν, μ] += -k*cft*PhiG*((Δp +Δm)*ixa(Kg, ν, μ)/Kg[3] - 2.0*eps_lut[μ, ν, 3]*Δm) - 1.0im*k*PhiR*AI*ixa(Nv, ν, μ)
            end
        end
    end
    return Ge, Gh
end


function Ewald2D(p::Int, Rv::Vector{T}, k::Union{T, Complex{T}}, kp::Vector{T}, lattice::lattice2D; Nmax::Int = 5) where  T<:Union{Float64, Int}
  Ra::Array{T} = lattice.R
  Ga::Array{T} = lattice.G
  AA::T = lattice.AWZC
  η::T = max(1.8*sqrt(pi/AA), 0.5*k/sqrt(6))
  η2::T = η*η
  Γp2 = gamma(p/2)
  cft::T = 0.25*AA/pi
  gv = similar(kp)
  Nv = similar(Rv)
  Σp  = 0.0
  N = 4*Nmax^2 + 1
  for n = 1:N
      inds = spiral(n)
      #Initiating reciprocal sum parameters
      for l = 1:2
          gv[l] = kp[l] + inds[1].*Ga[1,l] + inds[2].*Ga[2,l]
          Nv[l] = Rv[l] - inds[1].*Ra[1,l] - inds[2].*Ra[2,l]
      end
      g = norm(gv)
      r = norm(Nv)
      Σp += cis(dot(gv, Rv))*sf_gamma_inc(p/2, g*g*η2)/g^p + cft*cis(dot(kp, Nv))*sf_gamma_inc(1.0-p/2, 0.25*r*r/η2)*(0.5*r)^(p-2.0)
  end
  Σp = Σp/Γp2
  return Σp
end


#= Test the lattice sums

using PyPlot
pygui(true)
θ = pi/4
ϕ = 0
a = 1000
lat = lattice2D(a, a, pi/2)
nlam1 = 5000
λ1 = range(300, 5000, length = nlam1)
G1_E = zeros(ComplexF64, nlam1, 3, 3)
G1_H = zeros(ComplexF64, nlam1, 3, 3)
G2_E = zeros(ComplexF64, nlam1, 3, 3)
G2_H = zeros(ComplexF64, nlam1, 3, 3)

for (i, λi) = enumerate(λ1)
  k0 =  2*π/λi
  kp = [k0*sin(θ)*cos(ϕ), k0*sin(θ)*sin(ϕ)]
  k1 =  k0*sqrt(4)
  kp1 = [k1*sin(θ)*cos(ϕ), k1*sin(θ)*sin(ϕ)]
  G1_E[i, :, :], G1_H[i, :, :] = Lsum(zeros(3), k0, kp, lat, Nmax = 5)
  G2_E[i, :, :], G2_H[i, :, :] = Lsum(zeros(3), k1, kp1, lat, Nmax = 5)
end

begin
    plot(λ1, real(G1_H[:, 2, 3]), color = "Orange")
    plot(λ1, real(G2_H[:, 2, 3]), color = "black")
end
=#
