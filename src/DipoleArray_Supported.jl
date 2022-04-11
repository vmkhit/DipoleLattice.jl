
function zsqrt(x)
  #  ZSQRT - Square root for complex argument.
  #  Chose sign such that imaginary part is always positive.
  y = sqrt(x + 0.0im)
  y = y.*sign(imag(y + 1im*eps(Float64)))
end

function get_propagating_gvectors(k::Real, kx::Real, ky::Real, L::Lattice2D)
  Gmax =  round(Int, k/min(π/norm(mlat.R[1, :]), π/norm(mlat.R[2, :]))) + 1
  Gv = Number[]
  ng = 0
  for m = - Gmax:Gmax
    for n = - Gmax:Gmax
      gx = kx + m*L.G[1, 1] + n*L.G[2, 1]
      gy = ky + m*L.G[1, 2] + n*L.G[2, 2]
      if (gx*gx + gy*gy) < k*k
        append!(Gv, complex(gx, gy))
        ng = ng + 1
      end
    end
  end
  if ng < Gmax*Gmax
    nothing
  else
    throw("Gmax is  not sufficient for far field")
  end
  Gv = sort(Gv, by = x -> (abs(x), angle(x)), rev=false)  # put rev = true to go from large to small
  return ng, [real.(Gv) imag.(Gv)]
end

function Make_alpha3x3(α::Number; alpha3x3::Array{<:Number, 2} = fill(complex(1e-25, 1e-25), 3, 3))
  for i = 1:3
      alpha3x3[i, i] = α
  end
  return alpha3x3
end

function Make_alpha6x6(αee::Number, αmm::Number; alpha6x6::Array{<:Number, 2} = fill(complex(1e-25, 2e-15), 6, 6))
  for i = 1:3
      alpha6x6[i, i] = αee
      alpha6x6[i+3, i+3] = αmm
  end
  return alpha6x6
end

mutable struct Dipoles
  Centers::Array{<:Real}
  Alphas::AbstractVector{<:Array{<:Number}}
end

mutable struct Spheres
  Centers::Array{<:Real}
  Radii::Array{<:Real}
  EpsMu::Array{<:Number}
end


#======================  Utilities for S-Matrices =============================#

function Pair(NG::Integer, SA1::Array, SA2::Array, SA3::Array, SA4::Array,
                           SB1::Array, SB2::Array, SB3::Array, SB4::Array)
  II = Matrix(I, 2*NG, 2NG)
  A = inv(II - SA2*SB3)
  B = inv(II - SB3*SA2)
  SAB1 = SB1*A*SA1
  SAB2 = SB2 + SB1*SA2*B*SB4
  SAB3 = SA3 + SA4*SB3*A*SA1
  SAB4 = SA4*B*SB4
  return SAB1, SAB2, SAB3, SAB4
end


function SubstrateSMatrix(k::Real, qvec::Vector{<:Real}, NG::Integer, GVEC::Array{<:Real, 2}, S::Stack)
  S1 = zeros(Complex{Float64}, (2*NG, 2*NG))
  S2 = zeros(Complex{Float64}, (2*NG, 2*NG))
  S3 = zeros(Complex{Float64}, (2*NG, 2*NG))
  S4 = zeros(Complex{Float64}, (2*NG, 2*NG))
  R = zeros(Complex{Float64}, (2, 2))
  T = zeros(Complex{Float64}, (2, 2))
  S1 = revert_stack(S)
  for ii = 1:NG
    qG = norm()
    for i = 1:2
      T[1, i], R[1, i] = tmm_matrix(i, 2π/k, qvec + GVEC[:, ii], S)     # r, t from top
      T[2, i], R[2, i] = tmm_matrix(i, 2π/k, qvec + GVEC[:, ii], S1)    # r, t from bottom
    end
    for s = 1:2
      is = 2*(ii - 1) + s
      S1[is, is] = T[1, s]
      S2[is, is] = R[2, s]
      S3[is, is] = R[1, s]
      S4[is, is] = T[2, s]
    end
  end
  return S1, S2, S3, S4
end

function DipoleSmatrix(k::Real, qvec::Vector{<:Real}, NG::Integer, GVEC::Array{<:Real, 2}, D::Dipoles, L::Lattice2D)

  # constrcut the incidence wavector component
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  q  = sqrt(kx*kx + ky*ky)
  kz = k*cos(θ)
  cft = 2.0im*π*k1^2/(L.Area*kz)

  # define polarization unit - vectors
  es  = [ sin(ϕ),  cos(ϕ), 0.0]
  epm = [-cos(ϕ), -sin(ϕ), 0.0]
  epp = [ cos(ϕ),  sin(ϕ), 0.0]

  if θ > 0
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k1), -ky*kz/(q*k1), -q/k1]
    epp = [ kx*kz/(q*k1),  ky*kz/(q*k1), -q/k1]
  end
  uz = [0, 0, 1]           # unit vector along z - axis

  # incident field amplitudes
  F0s_p = [es; -epp];  F0s_m = [es; -epm]
  F0p_p = [epp; es];   F0p_m = [epm; es]

  # Allocate the memory
  M = zeros(ComplexF64, 6*Ns, 6*Ns)
  Fs_m = zeros(ComplexF64, 6*Ns)
  Fp_m = zeros(ComplexF64, 6*Ns)
  Fs_p = zeros(ComplexF64, 6*Ns)
  Fp_p = zeros(ComplexF64, 6*Ns)
  Ψm = zeros(ComplexF64, 6*Ns)
  Ψp = zeros(ComplexF64, 6*Ns)

  Gk = DyadSum(k1, zeros(3), [kx, ky], L, "GH")


  ν = 1
  for l = 1:6:6*Ns
    rν = D.Centers[ν, :]
    ψm = cis(kx*rν[1] + ky*rν[2] - kz*rν[3])
    ψp = cis(kx*rν[1] + ky*rν[2] + kz*rν[3])
    # Construct the external field at each dipole
    for s = 0:5
      Ψm[l+s] = ψm
      Ψp[l+s] = ψp
      Fs_m[l+s] = F0s_m[s+1]
      Fp_m[l+s] = F0p_m[s+1]
      Fs_p[l+s] = F0s_m[s+1]
      Fp_p[l+s] = F0p_m[s+1]
    end
    # construct the matrix M (see the Notes)
    invα = inv(D.Alphas[ν])
    Gk_ref = RDyadSum(2*π/k0, [kx, ky], rν, rν, L, S, 3)

    νp =  1
    for lp = 1:6:6*Ns
      rνp = D.Centers[νp, :]

      Gk_ννp = DyadGreenFn(k1, rν - rνp, [kx, ky], L, "GH")
      Gk_ννp_ref = RDyadSum(2*π/k0, [kx, ky], rνp, rν, L, S, 3)

      if ν == νp
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = invα[s+1, q+1] - Gk[s+1, q+1] #- Gk_ref[s+1, q+1]
        end
      else
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = - Gk_ννp[s+1, q+1] #- Gk_ννp_ref[s+1, q+1]
        end
      end
      νp = νp + 1
    end
    ν = ν + 1
  end

  Ω = inv(M)
  P_s_m = Ω*(Fs_m.*Ψm)
  P_p_m = Ω*(Fp_m.*Ψm)
  P_s_p = Ω*(Fs_p.*Ψp)
  P_p_p = Ω*(Fp_p.*Ψp)

end
