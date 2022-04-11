include("TMM.jl")
include("Lattice.jl")

function zsqrt(x)
  #  ZSQRT - Square root for complex argument.
  #  Chose sign such that imaginary part is always positive.
  y = sqrt(x + 0.0im)
  y = y.*sign(imag(y + 1im*eps(Float64)))
end

function RDyadSum(λ::Real, kp::Vector{<:Real}, r1::Vector{<:Real}, r2::Vector{<:Real}, L::Lattice2D, S::Stack, Kmax)
  n1 = S.Layers[1].mat.nk
  k = 2.0*pi*n1/λ

  Gee = zeros(ComplexF64, 3, 3)
  Gem = zeros(ComplexF64, 3, 3)
  Gme = zeros(ComplexF64, 3, 3)
  Gmm = zeros(ComplexF64, 3, 3)

  es =  [ 0.0, 1.0, 0.0]
  epp = [ 1.0, 0.0, 0.0]
  epm = [-1.0, 0.0, 0.0]

  for m = -Kmax:Kmax
    for n = -Kmax:Kmax
      gx = kp[1] + m*L.G[1, 1] +  n*L.G[2, 1]
      gy = kp[2] + m*L.G[1, 2] +  n*L.G[2, 2]
      q = sqrt(gx*gx + gy*gy)
      kzg = zsqrt(k*k - q*q)

      if q != 0
        es =  [-gy/q, gx/q, 0.0]
        epm = [-kzg*gx/(k*q), -kzg*gy/(k*q), -q/k]
        epp = [ kzg*gx/(k*q),  kzg*gy/(k*q), -q/k]
      end

      rs = tmm_matrix(1, λ, [gx, gy], S)[1]
      rp = tmm_matrix(2, λ, [gx, gy], S)[1]

      cft = 2.0im*pi*k^2/(kzg*L.Area)
      ψ = cis(gx*(r1[1] - r2[1]) + gy*(r1[2] - r2[2]) + kzg*abs(r1[3]+ r1[3]))

      rs = cft*rs*ψ
      rp = cft*rp*ψ

      Gee +=  rs.*kron(es,   es') + rp.*kron(epp, epm')
      Gem += -rs.*kron(es,  epm') + rp.*kron(epp,  es')
      Gme += -rs.*kron(epp,  es') + rp.*kron(es,  epm')
      Gmm +=  rs.*kron(epp, epm') + rp.*kron(es,   es')
    end
  end
  return [Gee Gem; Gme Gmm]
end



#=======================  Example of usage ====================================
using PyPlot
pygui(true)
include("utils.jl")
include("DyadSum.jl")


mat0 = Material()
mat1 = Material("epsmu",  2.1, 1.0)
mat2 = Material("epsmu",  2.1, 1.0)
L1 = Layer(mat0, 0)
L3 = Layer(mat1, 600)
L5 = Layer(mat1, 0)

a = 1000
mlat = lattice2D(a, a, pi/3)


θ = pi/3
nl = 6000
lam = range(300, 3000, length = nl)
RT = Array{Float64}(undef, nl, 2)

G_R = zeros(ComplexF64, nl, 6, 6)
G_0 = zeros(ComplexF64, nl, 6, 6)


for (l, λ) = enumerate(lam)
    k0 =  2π/λ
    kx = k0*sin(θ)

    ϵm = eps_Drude(λ, 4.0, 9.01, 0.1)
    mat3 = Material("epsmu", ϵm, 1.0)
    L2 = Layer(mat3, 20)
    L4 = Layer(mat3, 100)
    S = Stack([L1, L2, L3, L4, L5], zeros(4))

    RT[l, 1], RT[l, 2] = RT_calc(2, λ, [kx, 0.0], S)
    G_R[l, :, :] = RDyadSum(λ, [kx, 0], [0, 0, 100], [0, 0, 200], mlat, S, 20)
    G_0[l, :, :] = DyadSum(k0,  [0, 0, 100], [kx, 0], mlat, "GH")
end

plot(lam, RT[:, 1:2])

begin
    #plot(lam, RT[1, :], label = "R")
    #plot(lam, RT[2, :], label = "T")
    #plot(lam, 1 .- RT[1, :] - RT[2, :], label = "A")
    for i = 1:1
      for j = 1:2
        plot(lam, G_R[:, i, j] +  G_0[:, i, j], color = "black")
      end
    end
    legend()
end

=#
