"""
zqrt(z) returns square root of complex number x and corectly handles the Riemann branch
"""
function zsqrt(z)
  y = sqrt(z + 0.0im)
  y = y.*sign(imag(y + 1im*eps(Float64)))
end

function get_propagating_gvectors(k::Real, kx::Real, ky::Real, L::Lattice2D)
  Gmax =  round(Int, k/min(π/norm(mlat.R[1]), π/norm(mlat.R[2]))) + 1
  Gv = Number[]
  ng = 0
  for m = - Gmax:Gmax
    for n = - Gmax:Gmax
      gx = kx + m*L.G[1][1] + n*L.G[2][1]
      gy = ky + m*L.G[1][2] + n*L.G[2][2]
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

function getpolarizationvectors(ϵ, k0, θ, ϕ)
  # constrcut the incidence wavector components
  k1 = k0*sqrt(ϵ)
  kx = k0*sin(θ)*cos(ϕ)
  ky = k0*sin(θ)*sin(ϕ)
  q  = sqrt(kx*kx + ky*ky)
  kz = zsqrt(k1*k1 - q*q)
  zu = [0, 0, 1]

  if θ <= 1e-10   # care must be taken for normal incidence to avoid infinities
    Kp =  [0, 0,  1]
    Km =  [0, 0, -1]
    es  = [ sin(ϕ), cos(ϕ), 0.0]
    epm = [-cos(ϕ), sin(ϕ), 0.0]
    epp = [ cos(ϕ), sin(ϕ), 0.0]
  else
    Kp  = [kx, ky,  kz]
    Km  = [kx, ky, -kz]
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k1), -ky*kz/(q*k1), -q/k1]
    epp = [ kx*kz/(q*k1),  ky*kz/(q*k1), -q/k1]
    #es  = [-ky/q, kx/q, 0.0]
    #epm = [-kx*kz/(q*k1), -ky*kz/(q*k1), q/k1]
    #epp = [ kx*kz/(q*k1),  ky*kz/(q*k1), q/k1]
  end
  # We can check to see that they satisfy to
  #II1 = Kp * Kp' + epp * epp' + es*es'
  #II2 = Km * Km' + epm * epm' + es*es'
  #Amx = IxA(Km)
  #Bmx = (es* epm') - (epm*es')
  #Apx = IxA(Kp)
  #Bpx = (es* epp') - (epp*es')
  #abs.(Amx - Bmx) = 0
  #abs.(Apx - Bpx) = 0
  return Km, Kp, es, epm, epp
end

function RT_PM_MultiSpheres_v1(p_ol::String, k::Real, θ::Real, ϕ::Real, ϵμh::Vector{<:Real}, S::Spheres, L::Lattice2D)
  @assert pol == "s" || pol == "p"  # polarization should be "p" or "s" otherwise the function will throw error
  Ns = length(S.Radii)  # number of spheres in the unit cell

  # constrcut the incidence wavector components
  k  = k*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = 2.0im*π*k^2/(L.AreaWZC*kz)

  # define polarization unit - vectors and incidence field amplitudes
  if θ <= 1e-10   # care must be taken for normal incidence to avoid infinities
    es  = [ sin(ϕ), cos(ϕ), 0.0]
    epm = [-cos(ϕ), sin(ϕ), 0.0]
    epp = [ cos(ϕ), sin(ϕ), 0.0]
  else
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k), -ky*kz/(q*k), -q/k]
    epp = [ kx*kz/(q*k),  ky*kz/(q*k), -q/k]
  end

  # Electric fields for given polarization and propagation directtion
  Fs_p = [es; -epp];  Fs_m = [es; -epm]
  Fp_p = [epp; es];   Fp_m = [epm; es]

  # Allocate the memory
  M = zeros(ComplexF64, 6*Ns, 6*Ns)
  F = zeros(ComplexF64, 6*Ns)
  ν = 1
  for l = 1:6:6*Ns
    rν = S.Centers[ν, :]
    # External field phase at each dipole
    phiν = cis(kx*rν[1] + ky*rν[2] - kz*rν[3])
    # Fill in the external field vector
    for q = 0:5
      if pol == "s"
        F[l+q] = Fs_m[q+1]*phiν
      else
        F[l+q] = Fp_m[q+1]*phiν
      end
    end
    # construct the matrix M (see the Notes)
    ale, alm = al_em(2*π/k, [S.EpsMu[ν, 1], ϵμh[1]], [S.EpsMu[ν, 2], ϵμh[2]], [S.Radii[ν]])   # calculate sphere polarizabilities
    invα = Make_alpha6x6(1.0/ale, 1.0/alm)                                            # inverse of diagonal polarizability matrix
    νp =  1
    for lp = 1:6:6*Ns
      rνp = S.Centers[νp, :]
      Gk_ννp  =  DyadSum(k, rν - rνp, [kx, ky], L, "GH")
      GR_ννp  =  GHdyad(rν - rνp, k)
      if ν == νp
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = invα[s+1, q+1] - Gk_ννp[s+1, q+1]
        end
      else
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = -GR_ννp[s+1, q+1] - Gk_ννp[s+1, q+1]
        end
      end
      νp = νp + 1
    end
    ν = ν + 1
  end
  # calculate the induced dipoles
  P = inv(M)*F
  # allocate the matrix for total reflection coeficient r_{sνν'}, such that
  # if ν = "pol" = "s" => ν' = "p" and vice versa
  r_s = complex(0.0, 0.0)
  r_p = complex(0.0, 0.0)
  t_s = complex(0.0, 0.0)
  t_p = complex(0.0, 0.0)
  ν = 1
  for l = 1:6:6*Ns
    rν = S.Centers[ν, :]
    pm = P[l:l+5]

    ψν_p = cis(-kx*rν[1] - ky*rν[2] + kz*rν[3])
    ψν_m = cis(-kx*rν[1] - ky*rν[2] - kz*rν[3])

    r_s += dot(Fs_p, pm)*ψν_p
    r_p += dot(Fp_p, pm)*ψν_p
    t_s += dot(Fs_m, pm)*ψν_m
    t_p += dot(Fp_m, pm)*ψν_m
    ν = ν + 1
  end
  if pol == "s"
    t_s = 1 + cft*t_s
    t_p = cft*t_p
  else
    t_s = cft^t_s
    t_p = 1 + cft*t_p
  end
  r_s = cft*r_s; r_p = cft*r_p
  return r_s, r_p, t_s, t_p
end

function RT_PM_MultiSpheres_v2(pol::String, k::Real, θ::Real, ϕ::Real, ϵμh::Vector{<:Real}, S::Spheres, L::Lattice2D)
  @assert pol == "s" || pol == "p"  # polarization should be "p" or "s" otherwise the function will throw error
  Ns = length(S.Radii)  # number of spheres in the unit cell

  # constrcut the incidence wavector components
  k  = k*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = 2.0im*π*k^2/(L.AreaWZC*kz)
  # define polarization unit - vectors and incidence field amplitudes
  if θ <= 1e-10   # care must be taken for normal incidence to avoid infinities
    es  = [ sin(ϕ), cos(ϕ), 0.0]
    epm = [-cos(ϕ), sin(ϕ), 0.0]
    epp = [ cos(ϕ), sin(ϕ), 0.0]
  else
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k), -ky*kz/(q*k), -q/k]
    epp = [ kx*kz/(q*k),  ky*kz/(q*k), -q/k]
  end

  # Electric fields for given polarization and propagation directtion
  Fs_p = [es; -epp];  Fs_m = [es; -epm]
  Fp_p = [epp; es];   Fp_m = [epm; es]

  # Allocate the memory
  M = zeros(ComplexF64, 6*Ns, 6*Ns)
  F = zeros(ComplexF64, 6*Ns)
  Gk = DyadSum(k, zeros(3), [kx, ky], L, "GH")

  ν = 1
  for l = 1:6:6*Ns
    rν = S.Centers[ν, :]
    ψν = cis(kx*rν[1] + ky*rν[2] - kz*rν[3])
    # Construct the external field at each dipole
    for s = 0:5
      if pol == "s"
        F[l+s] = Fs_m[s+1].*ψν
      else
        F[l+s] = Fp_m[s+1].*ψν
      end
    end
    # construct the matrix M (see the Notes)
    ale, alm = al_em(2*π/k, [S.EpsMu[ν, 1], ϵμh[1]], [S.EpsMu[ν, 2], ϵμh[2]], [S.Radii[ν]])   # calculate sphere polarizabilities
    invα = Make_alpha6x6(1.0/ale, 1.0/alm)     # inverse of diagonal polarizability matrix
    νp =  1
    for lp = 1:6:6*Ns
      rνp = S.Centers[νp, :]
      Gk_ννp = DyadGreenFn(k, rν - rνp, [kx, ky], L, "GH")
      if ν == νp
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = invα[s+1, q+1] - Gk[s+1, q+1]
        end
      else
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = - Gk_ννp[s+1, q+1]
        end
      end
      νp = νp + 1
    end
    ν = ν + 1
  end
  # calculate the induced dipoles
  P = inv(M)*F
  # allocate the matrix for total reflection coeficient r_{sνν'}, such that
  # if ν = "pol" = "s" => ν' = "p" and vice versa

  r_s = complex(0.0, 0.0)
  r_p = complex(0.0, 0.0)
  t_s = complex(0.0, 0.0)
  t_p = complex(0.0, 0.0)
  ν = 1
  for l = 1:6:6*Ns
    rν = S.Centers[ν, :]
    pm = P[l:l+5]
    ψν_p = cis(-kx*rν[1] - ky*rν[2] + kz*rν[3])
    ψν_m = cis(-kx*rν[1] - ky*rν[2] - kz*rν[3])
    r_s += dot(Fs_p, pm)*ψν_p
    r_p += dot(Fp_p, pm)*ψν_p
    t_s += dot(Fs_m, pm)*ψν_m
    t_p += dot(Fp_m, pm)*ψν_m
    ν = ν + 1
  end
  if pol == "s"
    t_s = 1 + cft*t_s
    t_p = cft*t_p
  else
    t_s = cft*t_s
    t_p = 1 + cft*t_p
  end
  r_s = cft*r_s; r_p = cft*r_p
  return r_s, r_p, t_s, t_p
end

function RT_PM_DipoleCluster(pol::String, k::Real, θ::Real, ϕ::Real, ϵμh::Vector{<:Real}, S::Dipoles, L::Lattice2D)
  @assert pol == "s" || pol == "p"  # polarization should be "p" or "s" otherwise the function will throw error
  Ns = size(S.Centers, 1)  # number of spheres in the unit cell

  # constrcut the incidence wavector components
  k  = k*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = 2.0im*π*k^2/(L.AreaWZC*kz)

  # define polarization unit - vectors and incidence field amplitudes
  if θ <= 1e-10   # care must be taken for normal incidence to avoid infinities
    es  = [ sin(ϕ), cos(ϕ), 0.0]
    epm = [-cos(ϕ), -sin(ϕ), 0.0]
    epp = [ cos(ϕ), sin(ϕ), 0.0]
  else
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k), -ky*kz/(q*k), -q/k]
    epp = [ kx*kz/(q*k),  ky*kz/(q*k), -q/k]
  end

  # Electric fields for given polarization and propagation directtion
  Fs_p = [es; -epp];  Fs_m = [es; -epm]
  Fp_p = [epp; es];   Fp_m = [epm; es]

  # Allocate the memory
  M = zeros(ComplexF64, 6*Ns, 6*Ns)
  F = zeros(ComplexF64, 6*Ns)
  Gk = DyadSum(k, zeros(3), [kx, ky], L, "GH")
  ν = 1
  for l = 1:6:6*Ns
    rν = S.Centers[ν, :]
    ψν = cis(kx*rν[1] + ky*rν[2] - kz*rν[3])
    # Construct the external field at each dipole
    for s = 0:5
      if pol == "s"
        F[l+s] = Fs_m[s+1].*ψν
      else
        F[l+s] = Fp_m[s+1].*ψν
      end
    end
    # construct the matrix M (see the Notes)
    invα = inv(S.Alphas[ν])
    νp =  1
    for lp = 1:6:6*Ns
      rνp = S.Centers[νp, :]
      Gk_ννp = DyadGreenFn(k, rν - rνp, [kx, ky], L, "GH")
      if ν == νp
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = invα[s+1, q+1] - Gk[s+1, q+1]
        end
      else
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = - Gk_ννp[s+1, q+1]
        end
      end
      νp = νp + 1
    end
    ν = ν + 1
  end
  # calculate the induced dipoles
  P = inv(M)*F
  # allocate the matrix for total reflection coeficient r_{sνν'}, such that
  # if ν = "pol" = "s" => ν' = "p" and vice versa

  r_s = complex(0.0, 0.0)
  r_p = complex(0.0, 0.0)
  t_s = complex(0.0, 0.0)
  t_p = complex(0.0, 0.0)
  ν = 1
  for l = 1:6:6*Ns
    rν = S.Centers[ν, :]
    pm = P[l:l+5]
    ψν_p = cis(-kx*rν[1] - ky*rν[2] + kz*rν[3])
    ψν_m = cis(-kx*rν[1] - ky*rν[2] - kz*rν[3])
    r_s += dot(Fs_p, pm)*ψν_p
    r_p += dot(Fp_p, pm)*ψν_p
    t_s += dot(Fs_m, pm)*ψν_m
    t_p += dot(Fp_m, pm)*ψν_m
    ν = ν + 1
  end
  if pol == "s"
    t_s = 1 + cft*t_s
    t_p = cft*t_p
  else
    t_s = cft*t_s
    t_p = 1 + cft*t_p
  end
  r_s = cft*r_s; r_p = cft*r_p
  return r_s, r_p, t_s, t_p
end

function RT_P_DipoleCluster(pol::String, k::Real, θ::Real, ϕ::Real, ϵμh::Vector{<:Real}, S::Dipoles, L::Lattice2D)
  @assert pol == "s" || pol == "p"  # polarization should be "p" or "s" otherwise the function will throw error
  Ns = size(S.Centers, 1)  # number of spheres in the unit cell

  # constrcut the incidence wavector components
  k  = k*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = 2.0im*π*k^2/(L.AreaWZC*kz)

  # define polarization unit - vectors and incidence field amplitudes
  if θ <= 1e-10   # care must be taken for normal incidence to avoid infinities
    es  = [ sin(ϕ), cos(ϕ), 0.0]
    epm = [-cos(ϕ), sin(ϕ), 0.0]
    epp = [ cos(ϕ), sin(ϕ), 0.0]
  else
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k), -ky*kz/(q*k), -q/k]
    epp = [ kx*kz/(q*k),  ky*kz/(q*k), -q/k]
  end

  # Electric fields for given polarization and propagation directtion

  # Allocate the memory
  M = zeros(ComplexF64, 3*Ns, 3*Ns)
  F = zeros(ComplexF64, 3*Ns)
  Gk = DyadSum(k, zeros(3), [kx, ky], L, "GH")
  ν = 1
  for l = 1:3:3*Ns
    rν = S.Centers[ν, :]
    ψν = cis(kx*rν[1] + ky*rν[2] - kz*rν[3])
    # Construct the external field at each dipole
    for s = 0:2
      if pol == "s"
        F[l+s] = es[s+1].*ψν
      else
        F[l+s] = epm[s+1].*ψν
      end
    end
    # construct the matrix M (see the Notes)
    invα = inv(S.Alphas[ν])
    νp =  1
    for lp = 1:3:3*Ns
      rνp = S.Centers[νp, :]
      Gk_ννp = DyadGreenFn(k, rν - rνp, [kx, ky], L, "G")
      if ν == νp
        for s = 0:2, q = 0:2
          M[l+s, lp+q] = invα[s+1, q+1] - Gk[s+1, q+1]
        end
      else
        for s = 0:2, q = 0:2
          M[l+s, lp+q] = - Gk_ννp[s+1, q+1]
        end
      end
      νp = νp + 1
    end
    ν = ν + 1
  end
  # calculate the induced dipoles
  P = inv(M)*F
  # allocate the matrix for total reflection coeficient r_{sνν'}, such that
  # if ν = "pol" = "s" => ν' = "p" and vice versa

  r_s = complex(0.0, 0.0)
  r_p = complex(0.0, 0.0)
  t_s = complex(0.0, 0.0)
  t_p = complex(0.0, 0.0)
  ν = 1
  for l = 1:3:3*Ns
    rν = S.Centers[ν, :]
    p = P[l:l+2]
    ψν_p = cis(-kx*rν[1] - ky*rν[2] + kz*rν[3])
    ψν_m = cis(-kx*rν[1] - ky*rν[2] - kz*rν[3])
    r_s += dot(es,  p)*ψν_p
    r_p += dot(epp, p)*ψν_p
    t_s += dot(es,  p)*ψν_m
    t_p += dot(epm, p)*ψν_m
    ν = ν + 1
  end
  if pol == "s"
    t_s = 1 + cft*t_s
    t_p = cft*t_p
  else
    t_s = cft*t_s
    t_p = 1 + cft*t_p
  end
  r_s = cft*r_s; r_p = cft*r_p
  return r_s, r_p, t_s, t_p
end

function RT_M_DipoleCluster(pol::String, k::Real, θ::Real, ϕ::Real, ϵμh::Vector{<:Real}, S::Dipoles, L::Lattice2D)
  @assert pol == "s" || pol == "p"  # polarization should be "p" or "s" otherwise the function will throw error
  Ns = size(S.Centers, 1)  # number of spheres in the unit cell

  # constrcut the incidence wavector components
  k  = k*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = 2.0im*π*k^2/(L.AreaWZC*kz)

  # define polarization unit - vectors and incidence field amplitudes
  if θ <= 1e-10   # care must be taken for normal incidence to avoid infinities
    es  = [ sin(ϕ), cos(ϕ), 0.0]
    epm = [-cos(ϕ), sin(ϕ), 0.0]
    epp = [ cos(ϕ), sin(ϕ), 0.0]
  else
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k), -ky*kz/(q*k), -q/k]
    epp = [ kx*kz/(q*k),  ky*kz/(q*k), -q/k]
  end

  # Electric fields for given polarization and propagation directtion

  # Allocate the memory
  M = zeros(ComplexF64, 3*Ns, 3*Ns)
  F = zeros(ComplexF64, 3*Ns)
  Gk = DyadSum(k, zeros(3), [kx, ky], L, "GH")
  ν = 1
  for l = 1:3:3*Ns
    rν = S.Centers[ν, :]
    ψν = cis(kx*rν[1] + ky*rν[2] - kz*rν[3])
    # Construct the external magnetic field at each dipole
    for s = 0:2
      if pol == "s"
        F[l+s] = -epm[s+1].*ψν
      else
        F[l+s] = es[s+1].*ψν
       end
    end
    # construct the matrix M (see the Notes)
    invα = inv(S.Alphas[ν])
    νp =  1
    for lp = 1:3:3*Ns
      rνp = S.Centers[νp, :]
      Gk_ννp = DyadGreenFn(k, rν - rνp, [kx, ky], L, "G")
      if ν == νp
        for s = 0:2, q = 0:2
          M[l+s, lp+q] = invα[s+1, q+1] - Gk[s+1, q+1]
        end
      else
        for s = 0:2, q = 0:2
          M[l+s, lp+q] = - Gk_ννp[s+1, q+1]
        end
      end
      νp = νp + 1
    end
    ν = ν + 1
  end
  # calculate the induced dipoles
  P = inv(M)*F
  # allocate the matrix for total reflection coeficient r_{sνν'}, such that
  # if ν = "pol" = "s" => ν' = "p" and vice versa

  r_s = complex(0.0, 0.0)
  r_p = complex(0.0, 0.0)
  t_s = complex(0.0, 0.0)
  t_p = complex(0.0, 0.0)
  ν = 1
  for l = 1:3:3*Ns
    rν = S.Centers[ν, :]
    m = P[l:l+2]
    ψν_p = cis(-kx*rν[1] - ky*rν[2] + kz*rν[3])
    ψν_m = cis(-kx*rν[1] - ky*rν[2] - kz*rν[3])
    r_s += -dot(epp, m)*ψν_p
    r_p +=  dot(es,  m)*ψν_p
    t_s += -dot(epm, m)*ψν_m
    t_p +=  dot(es,  m)*ψν_m
    ν = ν + 1
  end
  if pol == "s"
    t_s = 1 + cft*t_s
    t_p = cft*t_p
  else
    t_s = cft*t_s
    t_p = 1 + cft*t_p
  end
  r_s = cft*r_s; r_p = cft*r_p
  return r_s, r_p, t_s, t_p
end

function EHFAR_PM_DipoleCluster(pol::String, k::Real, θ::Real, ϕ::Real, ϵμh::Vector{<:Real}, S::Dipoles, L::Lattice2D, x0::Real, y0::Real, z0::Real)
  @assert pol == "s" || pol == "p"  # polarization should be "p" or "s" otherwise the function will throw error
  Ns = size(S.Centers, 1)  # number of spheres in the unit cell

  # constrcut the incidence wavector components
  k  = k*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = 2.0im*π*k^2/(L.AreaWZC*kz)

  # define polarization unit - vectors and incidence field amplitudes
  es  = [ sin(ϕ), cos(ϕ), 0.0]
  epm = [-cos(ϕ), sin(ϕ), 0.0]
  epp = [ cos(ϕ), sin(ϕ), 0.0]
  uz = [0, 0, 1]           # unit vector along z - axis

  ψ = cis(kx*x0 + ky*y0 + kz*z0)
  if pol == "s"
    Einc = es.*ψ  ; Hinc = -epm.*ψ
  else
    Einc = epm.*ψ ; Hinc = es.*ψ
  end

  # for other than normal incidence
  if θ > 1e-5
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k), -ky*kz/(q*k), -q/k]
    epp = [ kx*kz/(q*k),  ky*kz/(q*k), -q/k]
  end

  # Electric fields for given polarization and propagation directtion
  Fs_p = [es; -epp];  Fs_m = [es; -epm]
  Fp_p = [epp; es];   Fp_m = [epm; es]

  # Allocate the memory
  M = zeros(ComplexF64, 6*Ns, 6*Ns)
  F = zeros(ComplexF64, 6*Ns)
  Gk = DyadSum(k, zeros(3), [kx, ky], L, "GH")
  ν = 1
  for l = 1:6:6*Ns
    rν = S.Centers[ν, :]
    ψν = cis(kx*rν[1] + ky*rν[2] - kz*rν[3])
    # Construct the external field at each dipole
    for s = 0:5
      if pol == "s"
        F[l+s] = Fs_m[s+1].*ψν
      else
        F[l+s] = Fp_m[s+1].*ψν
      end
    end
    # construct the matrix M (see the Notes)
    invα = inv(S.Alphas[ν])
    νp =  1
    for lp = 1:6:6*Ns
      rνp = S.Centers[νp, :]
      Gk_ννp = DyadGreenFn(k, rν - rνp, [kx, ky], L, "GH")
      if ν == νp
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = invα[s+1, q+1] - Gk[s+1, q+1]
        end
      else
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = - Gk_ννp[s+1, q+1]
        end
      end
      νp = νp + 1
    end
    ν = ν + 1
  end
  # calculate the induced dipoles
  P = inv(M)*F
  # Now calculate the scattered fields including many orders at r = (0, 0, z0) plane

  Escat = zeros(ComplexF64, 2, 3)
  Hscat = zeros(ComplexF64, 2, 3)
  Kgm =  zeros(ComplexF64, 3)
  Kgp =  zeros(ComplexF64, 3)
  ng, Gv = get_propagating_gvectors(k, kx, ky, L)    # kp + gv lattice vectors that are incide the light cone

  for g = 1:ng
    kg_x,  kg_y = Gv[g, :]
    q = sqrt(kg_x*kg_x + kg_y*kg_y)
    kg_z = zsqrt(k*k - kg_x*kg_x - kg_y*kg_y)
    cftg = 2.0im*π/(kg_z*L.Area)
    Kgm[1] = kg_x;   Kgp[1] = kg_x
    Kgm[2] = kg_y;   Kgp[2] = kg_y
    Kgm[3] = -kg_z;  Kgp[3] = kg_z

    ν = 1
    for l = 1:6:6*Ns
      rν = S.Centers[ν, :]
      pν = P[l:l+2]
      mν = P[l+3:l+5]

      ψm = cis(kg_x*(x0 - rν[1]) + kg_y*(y0-rν[2]) + kg_z*abs(-z0 - rν[3]))
      ψp = cis(kg_x*(x0 - rν[1]) + kg_y*(y0-rν[2]) + kg_z*abs(z0 - rν[3]))

      Kp_dot_p = (kg_x*pν[1] + kg_y*pν[2] + kg_z*pν[3])
      Km_dot_p = (kg_x*pν[1] + kg_y*pν[2] - kg_z*pν[3])

      Kp_dot_m = (kg_x*mν[1] + kg_y*mν[2] + kg_z*mν[3])
      Km_dot_m = (kg_x*mν[1] + kg_y*mν[2] - kg_z*mν[3])

      for l = 1:3
        l1 = mod(l,   3) + 1
        l2 = mod(l+1, 3) + 1
        Escat[1, l] += cftg*(k*k*pν[l] - Kp_dot_p*Kgp[l] - k*(Kgp[l1]*mν[l2] - Kgp[l2]*mν[l1]))*ψp
        Escat[2, l] += cftg*(k*k*pν[l] - Km_dot_p*Kgm[l] - k*(Kgm[l1]*mν[l2] - Kgm[l2]*mν[l1]))*ψm
        Hscat[1, l] += cftg*(k*k*mν[l] - Kp_dot_m*Kgp[l] + k*(Kgp[l1]*pν[l2] - Kgp[l2]*pν[l1]))*ψp
        Hscat[2, l] += cftg*(k*k*mν[l] - Km_dot_m*Kgm[l] + k*(Kgm[l1]*pν[l2] - Kgm[l2]*pν[l1]))*ψm
      end
      ν = ν + 1
    end
  end
  return Einc, Hinc, Escat, Hscat
end


function RT_PM_DipoleCluster_diffractive(pol::String, k0::Real, θ::Real, ϕ::Real, ϵμh::Vector{<:Real}, S::Dipoles, L::Lattice2D)
  @assert pol == "s" || pol == "p"  # polarization should be "p" or "s" otherwise the function will throw error
  Ns = size(S.Centers, 1)  # number of spheres in the unit cell

  # constrcut the incidence wavector component
  nϵ1 = sqrt(ϵμh[1]*ϵμh[2])
  k1  = k0*nϵ1     # wavenumber of the medium

  kx = k1*sin(θ)*cos(ϕ)
  ky = k1*sin(θ)*sin(ϕ)
  q  = k1*sin(θ)
  kz = k1*cos(θ)

  cft = 2.0im*π*k1^2/(L.AreaWZC*kz)

  # define polarization unit - vectors and incidence field amplitudes
  es  = [ sin(ϕ), cos(ϕ), 0.0]
  epm = [-cos(ϕ), -sin(ϕ), 0.0]
  epp = [ cos(ϕ), sin(ϕ), 0.0]
  uz = [0, 0, 1]           # unit vector along z - axis

  if θ > 0
    es  = [-ky/q, kx/q, 0.0]
    epm = [-kx*kz/(q*k1), -ky*kz/(q*k1), -q/k1]
    epp = [ kx*kz/(q*k1),  ky*kz/(q*k1), -q/k1]
  end


  # Electric fields for given polarization and propagation directtion
  Fs_p = [es; -epp];  Fs_m = [es; -epm]
  Fp_p = [epp; es];   Fp_m = [epm; es]

  # Allocate the memory
  M = zeros(ComplexF64, 6*Ns, 6*Ns)
  F = zeros(ComplexF64, 6*Ns)
  Gk = DyadSum(k1, zeros(3), [kx, ky], L, "GH")

  ν = 1
  for l = 1:6:6*Ns
    rν = S.Centers[ν, :]
    ψν = cis(kx*rν[1] + ky*rν[2] - kz*rν[3])
    # Construct the external field at each dipole
    for s = 0:5
      if pol == "s"
        F[l+s] = Fs_m[s+1].*ψν
      else
        F[l+s] = Fp_m[s+1].*ψν
      end
    end
    # construct the matrix M (see the Notes)
    invα = inv(S.Alphas[ν])
    νp =  1
    for lp = 1:6:6*Ns
      rνp = S.Centers[νp, :]
      Gk_ννp = DyadGreenFn(k1, rν - rνp, [kx, ky], L, "GH")
      if ν == νp
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = invα[s+1, q+1] - Gk[s+1, q+1]
        end
      else
        for s = 0:5, q = 0:5
          M[l+s, lp+q] = - Gk_ννp[s+1, q+1]
        end
      end
      νp = νp + 1
    end
    ν = ν + 1
  end

  # calculate the induced dipoles
  P = inv(M)*F
  # Now calculate the scattered fields including many orders at r = (0, 0, z0) plane
  Kgm =  zeros(ComplexF64, 3)
  Kgp =  zeros(ComplexF64, 3)
  ng, Gv = get_propagating_gvectors(k1, kx, ky, L)    # kp + gv lattice vectors that are incide the light cone

  z0 = 20*π/k1  # large enough to consider it to be a far field  z0 = 10 λ

  ψ = cis(kz*z0)
  if pol == "s"
    Einc = es.*ψ  ; Hinc = -epm.*ψ
  else
    Einc = epm.*ψ ; Hinc = es.*ψ
  end
  Iinc = real(dot(cross(Einc, conj(Hinc)), -uz))   # incidence intensity

  R = T = 0.0
  for g = 1:ng
    kg_x,  kg_y = Gv[g, :]
    kg_z = zsqrt(k1*k1 - kg_x*kg_x - kg_y*kg_y)
    cftg = 2.0im*π/(kg_z*L.Area)

    Kgm[1] =  kg_x;   Kgp[1] = kg_x
    Kgm[2] =  kg_y;   Kgp[2] = kg_y
    Kgm[3] = -kg_z;   Kgp[3] = kg_z

    Escat = zeros(ComplexF64, 2, 3)
    Hscat = zeros(ComplexF64, 2, 3)

    if (kg_x == kx) && (kg_y == ky)
      Escat[2, :] += Einc
      Hscat[2, :] += Hinc
    end

    ν = 1
    for l = 1:6:6*Ns
      rν = S.Centers[ν, :]
      pν = P[l:l+2]
      mν = P[l+3:l+5]

      ψm = cis(-kg_x*rν[1] - kg_y*rν[2] + kg_z*abs(-z0 - rν[3]))
      ψp = cis(-kg_x*rν[1] - kg_y*rν[2] + kg_z*abs(z0 - rν[3]))

      Kp_dot_p = (kg_x*pν[1] + kg_y*pν[2] + kg_z*pν[3])
      Km_dot_p = (kg_x*pν[1] + kg_y*pν[2] - kg_z*pν[3])

      Kp_dot_m = (kg_x*mν[1] + kg_y*mν[2] + kg_z*mν[3])
      Km_dot_m = (kg_x*mν[1] + kg_y*mν[2] - kg_z*mν[3])

      for l = 1:3
        l1 = mod(l,   3) + 1
        l2 = mod(l+1, 3) + 1
        Escat[1, l] += cftg*(k1*k1*pν[l] - Kp_dot_p*Kgp[l] - k1*(Kgp[l1]*mν[l2] - Kgp[l2]*mν[l1]))*ψp
        Escat[2, l] += cftg*(k1*k1*pν[l] - Km_dot_p*Kgm[l] - k1*(Kgm[l1]*mν[l2] - Kgm[l2]*mν[l1]))*ψm
        Hscat[1, l] += cftg*(k1*k1*mν[l] - Kp_dot_m*Kgp[l] + k1*(Kgp[l1]*pν[l2] - Kgp[l2]*pν[l1]))*ψp
        Hscat[2, l] += cftg*(k1*k1*mν[l] - Km_dot_m*Kgm[l] + k1*(Kgm[l1]*pν[l2] - Kgm[l2]*pν[l1]))*ψm
      end
      # Onece the scattered fields are know we can constract total Reflectance and Transmittance from Pointing theorem
      ν = ν + 1
    end
    Rmn = real(dot(cross(Escat[1, :], conj(Hscat[1, :])),  uz))/Iinc
    Tmn = real(dot(cross(Escat[2, :], conj(Hscat[2, :])), -uz))/Iinc
    T += Tmn
    R += Rmn
  end
  return R, T, 1 - R - T
end
