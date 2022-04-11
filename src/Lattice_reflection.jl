####### ELECTRIC DIPOLE ONLY ##############################
"""
    rt_P_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αe::Array{<:Number, 2})
  This function returns the field rflection and transmission coeficienets of a 2D array of dipoles with only electric dipole moment.
  # Arguments
  - `pol::String`: Polarization of incidence plane wave
  - `λ::Real` : The free space wavelength of the incidence light
  - `θ::Real` : The azimuthal angle of the incidence plane wave
  - `ϕ::Real` : The polar angle of the the incidence plane wave
  - `ϵμh::Vector{<:Real}` : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].
  -  `lattice::Lattice2D` : Lattice objects that defines teh properties of 2D-lattice
  - `αe` : 3x3 electric polarizability tensor.
"""
function rt_P_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αe::Array{<:Number, 2})
    @assert pol == "s" || pol == "p"

    k0 = 2*π/λ
    # constrcut the incidence wavector components
    k  = k0*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
    kx = k*sin(θ)*cos(ϕ)
    ky = k*sin(θ)*sin(ϕ)
    kz = k*cos(θ)
    q  = k*sin(θ)
    cft = (2.0im*pi*k^2)/(kz*lattice.AreaWZC)

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

    if pol == "s"
      E0 = es
    else
      E0 = epm
    end
    # calculate the induced dipole moments
    Ge =   DyadSum(k, zeros(3), [kx, ky], lattice, "G", Nmax = 10)
    p  = inv(inv(αe) - Ge)*E0
    # calculate the r and t coeficients
    rs =  cft*dot(p, es)
    rp =  cft*dot(p, epp)
    ts = 1 + rs
    tp = 1 + cft*dot(p, epm)
    return rs, rp, ts, tp
end

#====================   Magnetic dipole ONLY ==================================#
"""
    rt_M_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αm::Array{<:Number, 2})
  This function returns the field rflection and transmission coeficienets of a 2D array of dipoles with only magnetic dipole moments.
  # Arguments
  - `pol::String`: Polarization of incidence plane wave
  - `λ::Real` : The free space wavelength of the incidence light
  - `θ::Real` : The azimuthal angle of the incidence plane wave
  - `ϕ::Real` : The polar angle of the the incidence plane wave
  - `ϵμh::Vector{<:Real}` : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].
  -  `lattice::Lattice2D` : Lattice objects that defines teh properties of 2D-lattice
  - `αm` : 3x3 magnetic polarizability tensor.
"""
function rt_M_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αm::Array{<:Number, 2})
  @assert pol == "s" || pol == "p"

  k0 = 2*π/λ
  # constrcut the incidence wavector components
  k  = k0*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = (2.0im*pi*k^2)/(kz*lattice.AreaWZC)
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

  if pol == "s"
    H0 = -epm    # note that this is different as compare to p dipole
  else
    H0 = es
  end

  Ge =   DyadSum(k, zeros(3), [ks, ky], lattice, "G", Nmax = 10)
  m = inv(inv(αm) - Ge)*H0

  rs =  cft*dot(-epp, m)
  rp =  cft*dot(es, m)
  ts =  1 + cft*dot(-epm, m)
  tp =  1 + rp
  return rs, rp, ts, tp
end

#================  Both ELECTRIC and MAGNETCI DIPOLES =========================#
"""
    rt_PM_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::lattice2D, αem::Array{<:Number, 2})
  This function returns the field rflection and transmission coeficienets of a 2D array of dipoles with both electric and magnetic dipole moments.
  # Arguments
  - `pol::String`: Polarization of incidence plane wave
  - `λ::Real` : The free space wavelength of the incidence light
  - `θ::Real` : The azimuthal angle of the incidence plane wave
  - `ϕ::Real` : The polar angle of the the incidence plane wave
  - `ϵμh::Vector{<:Real}` : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].
  -  `lattice::Lattice2D` : Lattice objects that defines teh properties of 2D-lattice
  - `αem` : 6x6 magneto-electrical polarizability tensor.
"""
function rt_PM_dipole(pol::String, λ::Real, θ::Real, ϕ::Real,  ɛμh::Vector{<:Real}, lattice::Lattice2D, α::Array{<:Number})
  @assert pol == "s" || pol == "p"
  k0 = 2*π/λ
  # constrcut the incidence wavector components
  k  = k0*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = (2.0im*pi*k^2)/(kz*lattice.AreaWZC)
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

  GEM = DyadSum(k, zeros(3), [kx, ky], lattice, "GH", Nmax = 10)

  if pol == "s"
    pm = inv(inv(α) - GEM)*Fs_m
    ts = 1 + cft.*dot(Fs_m, pm)
    tp = cft*dot(Fp_m, pm)
  else
    pm = inv(inv(α) - GEM)*Fp_m
    ts = cft.*dot(Fs_m, pm)
    tp = 1 + cft*dot(Fp_m, pm)
  end
  rs = cft*dot(Fs_p, pm)
  rp = cft*dot(Fp_p, pm)
  return rs, rp, ts, tp
end

#================  Both ELECTRIC and MAGNETCI DIPOLES Full r/t matrix =========#
"""
    rt_PM_dipole_FullMatrix(λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::lattice2D, α::Array{<:Number, 2})
  This function returns the field rflection and transmission coeficienets of a 2D array of dipoles with only electrica dipole moment.
  It is similar to `rt_PM_dipole` but returns the coeficient for both polarizations at once.
  # Arguments
  - `pol::String`: Polarization of incidence plane wave
  - `λ::Real` : The free space wavelength of the incidence light
  - `θ::Real` : The azimuthal angle of the incidence plane wave
  - `ϕ::Real` : The polar angle of the the incidence plane wave
  - `ϵμh::Vector{<:Real}` : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].
  -  `lattice::Lattice2D` : Lattice objects that defines teh properties of 2D-lattice
  - `α` : 6x6 electrical polarizability tensor.
"""
function rt_PM_dipole_FullMatrix(λ::Real, θ::Real, ϕ::Real,  ɛμh::Vector{<:Real}, lattice::Lattice2D, α::Array{<:Number})
  @assert pol == "s" || pol == "p"
  k0 = 2*π/λ
  # constrcut the incidence wavector components
  k  = k0*sqrt(ϵμh[1]*ϵμh[2])      # wavenumber of the medium
  kx = k*sin(θ)*cos(ϕ)
  ky = k*sin(θ)*sin(ϕ)
  kz = k*cos(θ)
  q  = k*sin(θ)
  cft = (2.0im*pi*k^2)/(kz*lattice.AreaWZC)

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

  GEM = DyadSum(k, zeros(3), [kx, ky], lattice, "GH", Nmax = 10)
  pms = inv(inv(α) - GEM)*Fs_m
  tss = 1 + cft*dot(Fs_m, pm_s)                    # incidence "s" - transmits "s"
  tsp = cft*dot(Fp_m, pm_s)                    # incidence "s" - transmits "p"
  rss = cft*dot(Fs_p, pm_s)                    # incidence "s" - reflects  "s"
  rsp = cft*dot(Fp_p, pm_s)                    # incidence "s" - reflects  "p"

  pm_p = inv(inv(α) - GEM)*Fp_m
  tps = cft*dot(Fs_m, pm_p)                        # incidence "p" - transmits "s"
  tpp = 1 + cft*dot(Fp_m, pm_p)                # incidence "p" - transmits "p"
  rps = cft*dot(Fs_p, pm_p)                    # incidence "p" - reflects  "s"
  rpp = cft*dot(Fp_p, pm_p)                    # incidence "p" - reflects  "p"

  return rss, rsp, rps, rpp, tss, tsp, tps, tpp
end
