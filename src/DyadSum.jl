using  SpecialFunctions
using LinearAlgebra

const eps_lut = cat([0 0 0; 0 0 1; 0 -1 0], [0 0 -1; 0 0 0; 1 0 0], [0 1 0; -1 0 0; 0 0 0], dims=3);
function IxA(a::Vector{<:Real}) where T<:Number
    M = zeros(eltype(a), (3,3))
    for j = 1:3
        for i = 1:3
            for k = 1:3
                M[i,j] += eps_lut[j, i, k]*a[k]
            end
        end
    end
    return M
end

function zsqrt(x)
  #  ZSQRT - Square root for complex argument.
  #  Chose sign such that imaginary part is always positive.
  y = sqrt(x + 0.0im)
  y = y.*sign(imag(y + 1im*eps(Float64)))
end

function βz(k::Number, q::Real)
  return zsqrt(k^2-q^2)
end


function g0(r::Real, k::Number)
  return cis(k*r)/r
end

function GL(r::Real, k::Number)
  return k^2 + im*k/r-one(r)/r^2
end
function GT(r::Real, k::Number)
  return 3.0/r^4-3.0im*k/r^3-k^2/r^2
end

function GxI(r::Real, k::Number)
  return im*k/r-one(r)/r^2
end

function A(r::Real, k::Number, η::Real)
  return 0.5*g0(r,k)*erfc(r*η+0.5im*k/η)
end

function B(r::Real, k::Number, η::Real)
  return 2.0*η*exp(-(r*η)^2+0.25*(k/η)^2)/sqrt(pi)
end

function C(z::Real, kz::Number, η::Real)
  return 0.5*cis(kz*z)*erfc(-z*η-0.5im*kz/η)/kz
end


# Electric Dyadic GF
function G0(rv::Vector{<:Real}, k::Number)
  r = norm(rv)
  return GL(r, k)*g0(r,k) .*diagm(0=>[1, 1, 1]) + GT(r,k)*g0(r,k).*kron(rv, rv')
end
#magnetic Dyadic GF
function H0(rv::Vector{<:Real}, k::Number)
  r = norm(rv)
  return im*k*GxI(r, k)*g0(r, k).*IxA(rv)
end

# 6x6 Dyadic Function
function GHdyad(rv::Vector{<:Real}, k::Number)
  G = G0(rv, k)
  H = H0(rv, k)
  return vcat(hcat(G, H), hcat(-H, G))
end

function S_singular(r::Real, k::Number, η::Real)
  if r == zero(r)
    return 2.0*η*exp(0.25*(k/η)^2)/sqrt(pi)+im*k*erfc(-0.5im*k/η)
  else
    return A(r, k, -η) - A(r,-k, η)
  end
end

function G_singular(rv::Vector{<:Real}, k::Number, η::Real)
  r = norm(rv)
  if r == 0.0
    return (4.0*η*(k^2-η^2)*exp(0.25*(k/η)^2)/(3.0*sqrt(pi)) + 2.0im*k^3*erfc(-0.5im*k/η)/3.0).*diagm(0=>[1, 1, 1])
  else
    return (GL(r,k)*A(r, k,-η)-GL(r,-k)*A(r,-k, η)+B(r, k, η)/r^2).*diagm(0=>[1, 1, 1]) -
           ((3.0/r^4+2.0*(η/r)^2)*B(r,k,η)-GT(r, k)*A(r,k,-η)+GT(r,-k)*A(r,-k,η)).*kron(rv, rv')
  end
end

function H_singular(rv::Vector{<:Real}, k::Number, η::Real)
  r = norm(rv)
  if r == zero(r)
    return zeros(eltype(r), 3, 3)
  else
    return im*k*(B(r, k, η)/r^2 + GxI(r,k)*A(r,k,-η) - GxI(r,-k)*A(r,-k, η)).*IxA(rv)
  end
end

function S_dir(rv::Vector{<:Real}, k::Number, kp::Vector{<:Real}, Rmn::Vector{<:Real}, η::Real)
  r = norm(rv-[Rmn; 0.0])
  return (A(r, k, η) + A(r,-k,η))*cis(kp'*Rmn)
end

function S_rec(rv::Vector{<:Real}, k::Number, kp::Vector{<:Real}, Gmn::Vector{<:Real}, η::Real)
  z = rv[3]
  Gv = [kp+Gmn; 0.0]
  kz = βz(k, norm(Gv))
  return (C(abs(z), kz, η) + C(-abs(z), kz, η))*cis(rv'*Gv)
end

function G_dir(rv::Vector{<:Real}, k::Number, kp::Vector{<:Real}, Rmn::Vector{<:Real}, η::Real)
  R = rv - [Rmn; 0.0]
  r = norm(R)
  return cis(dot(kp', Rmn)).*((GL(r,k)*A(r, k , η)+ GL(r, -k)*A(r, -k, η)- B(r, k, η)/r^2).*diagm(0=>[1, 1, 1])+
         (GT(r, k)*A(r, k, η) + GT(r, -k)*A(r, -k, η) + (3.0/r^4+2.0*(η/r)^2)*B(r, k, η)).*kron(R, R'))
end

function G_rec(rv::Vector{<:Real}, k::Number, kp::Vector{<:Real}, Gmn::Vector{<:Real}, η::Real)
  zu = [0.0, 0.0, 1.0]
  z = rv[3]
  Gv = [kp + Gmn; 0.0]
  kz = βz(k, norm(Gv))
  M1 = k*k.*diagm(0=>[1, 1, 1])-kron(Gv, Gv')-kz*kz.*kron(zu, zu') - kz.*(kron(Gv, zu') + kron(zu, Gv'))
  M2 = k*k.*diagm(0=>[1, 1, 1])-kron(Gv, Gv')-kz*kz.*kron(zu, zu') + kz.*(kron(Gv, zu') + kron(zu, Gv'))

  return (im*B(abs(z), kz, η).*kron(zu, zu') + C(abs(z), kz, η).*M1+ C(-abs(z), kz, η).*M2).*cis(dot(rv', Gv))
end

function H_dir(rv::Vector{<:Real}, k::Number, kp::Vector{<:Real}, Rmn::Vector{<:Real}, η::Real)
  R::Vector = rv - [Rmn; 0.0]
  r = norm(R)
  return -im*k*cis(kp'*Rmn)*(B(r, k, η)/r^2 - GxI(r, -k)*A(r, -k, η) - GxI(r, k)*A(r, k, η)).*IxA(R)
end

function H_rec(rv::Vector{<:Real}, k::Number, kp::Vector{<:Real}, Gmn::Vector{<:Real}, η::Real)
  zu = [0.0, 0.0, 1.0]
  z = rv[3]
  Gv = [kp + Gmn; 0.0]
  kz = βz(k, norm(Gv))
  return -k*cis(rv'*Gv).*((C(abs(z), kz, η)+C(-abs(z), kz, η)).*IxA(Gv)+kz*(C(abs(z), kz, η) - C(-abs(z), kz, η)).*IxA(zu))
end

function ScalarSum(k::Number, rv::Vector{<:Real}, kp::Vector{<:Real}, L::Lattice2D, Nmax::Integer = 5)
  η = 1.8*sqrt(pi/L.AreaWZC)
  cft = 2.0im*pi/L.AreaWZC
  R = L.R
  G = L.G
  SS = zero(ComplexF64)
  for m = -Nmax:Nmax
    for n = -Nmax:Nmax
      if n^2+m^2 != 0
        SS += S_dir(rv, k, kp, m.*R[1, :]+n.*R[2, :], η)
      end
      SS += cft*S_rec(rv, k, kp, m.*G[1, :]+n.*G[2, :], η)
    end
  end
  return SS - S_singular(norm(rv), k, η)
end

"""
  DyadSum(k::Number, rv::Vector{<:Real}, kp::Vector{<:Real}, L::Lattice2D, GorH::String; Nmax::Integer = 5)
"""
function DyadSum(k::Number, rv::Vector{<:Real}, kp::Vector{<:Real}, L::Lattice2D, GorH::String; Nmax::Integer = 5 )

  @assert GorH == "G" || GorH == "H" || GorH =="GH"

  η = 1.8*sqrt(pi/L.AreaWZC)
  η = max(η, k/2/sqrt(6))
  cft = 2.0im*pi/L.AreaWZC

  R = L.R
  G = L.G
  if GorH == "GH"
    SSH = zeros(ComplexF64, (3,3))
    SSG = zeros(ComplexF64, (3,3))
    for m = -Nmax:Nmax
      for n = -Nmax:Nmax
        if n^2 + m^2 != 0
          SSG +=  G_dir(rv, k, kp, m.*R[1, :]+n.*R[2, :], η)
          SSH += H_dir(rv, k, kp, m.*R[1, :]+n.*R[2, :], η)
        end
        SSG += cft.*G_rec(rv, k, kp, m.*G[1, :]+n.*G[2, :], η)
        SSH += cft.*H_rec(rv, k, kp, m.*G[1, :]+n.*G[2, :], η)
      end
    end
    SSG = SSG - G_singular(rv, k, η)
    SSH = SSH - H_singular(rv, k, η)
    return [SSG SSH; -SSH SSG]
  elseif GorH == "G"
    SS = zeros(ComplexF64, 3, 3)
    for m = -Nmax:Nmax
      for n = -Nmax:Nmax
        if n^2+m^2 !=0
           SS += G_dir(rv, k, kp, m.*R[1, :]+n.*R[2, :], η)
        end
        SS += cft.*G_rec(rv, k, kp, m.*G[1, :]+n.*G[2, :], η)
      end
    end
    SS = SS - G_singular(rv, k, η)
    return SS
  else
    SS = zeros(Complex{AbstractFloat}, (3, 3))
    for m = -Nmax:Nmax
      for n = -Nmax:Nmax
        if n^2+m^2 !=0
          SS += H_dir(rv, k, kp, m.*R[1, :]+n.*R[2, :], η)
        end
        SS += cft.*H_rec(rv, k, kp, m.*G[1, :]+n.*G[2, :], η)
      end
    end
    SS = SS - H_singular(rv, k, η)
    return SS
  end
end

function ScalarGreenFn(k::Number, rv::Vector{<:Real}, kp::Vector{<:Real}, L::Lattice2D, Nmax::Integer = 3)
  return ScalarSum(λ, ɛμ, rv, kp, L, Nmax) + g0(norm(rv), 2.0*pi*sqrt(prod(ɛμ))/λ)
end

function DyadGreenFn(k::Number, rv::Vector{<:Real}, kp::Vector{<:Real}, L::Lattice2D, GorH::String, Nmax::Integer = 3)
  @assert GorH == "G" || GorH == "H" || GorH == "GH"
  if GorH == "GH"
    return DyadSum(k, rv, kp, L, GorH, Nmax = Nmax) + [G0(rv, k) H0(rv, k); -H0(rv, k) G0(rv, k)]
  elseif GorH == "G"
    return DyadSum(k, rv, kp, L, GorH, Nmax = Nmax) + G0(rv, k)
  else
    return DyadSum(k, rv, kp, L, GorH, Nmax = Nmax) + H0(rv, k)
  end
end


#=Test the lattice sums
include("Lattice.jl")
lat = lattice2D(1.0, 1.0, pi/3)
k = 2*pi/600
kp = k.*[pi/4, 0]
rv = zeros(3)
η = 1.8*sqrt(pi/lat.Area)

S = DyadSum(k, rv, kp, lat, "H")
DyadGreenFn(k, rand(3), kp, lat, "G")
=#
