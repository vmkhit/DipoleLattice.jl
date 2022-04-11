using LinearAlgebra

"""
    epston(epsilon::Number)
This is a helper function to convert dielectric function to refractive index.
# Arguments
- `epsilon::Number`: Complex dielectric function.
# Example
```julia-repl
julia> n = epston(12 + 3im)
3.490653010315762 + 0.42971902264909145im
```
"""
function epston(epsilon::Number)
    n = sqrt(abs(epsilon) + real(epsilon))/sqrt(2)
    k = sqrt(abs(epsilon) - real(epsilon))/sqrt(2)
    return complex(n, k)
end


"""
    Material
This is struct object that defines material datastruct
# Arguments
- `id::String`: id="nk" or id = "epsmu" defines the input type. If the id="nk", we need to specify the refractive index. If id = "epsmu", the material will be defined through dielectric function and magnetic permeability.
- `eps::Number`: Dielectric function of the material
- `mu::Number`: magnetic permeability of the material
- `nk::Number`: complex refractive index of the material
# Examples
```julia-repl
julia> M1 = Material("nk", 2, 3)
Material("nk", -5.0 + 12.0im, 1.0, 2.0 + 3.0im)
julia> M2 = Material("epsmu", 2+3im, 1.0 + 0.0im)
Material("epsmu", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im)
```
"""
mutable struct Material
    id::String
    eps::Number
    mu::Number
    nk::Number
    function Material(id::String = "nk", kwargs...)
        @assert id == "epsmu" || id == "nk"
        if isempty(kwargs)
            eps = mu = nk = 1.0
        else
            if id == "nk"
                n = kwargs[1]
                k = kwargs[2]
                nk = n + 1.0im*k
                eps = nk^2
                mu = 1.0
            else
                if length(kwargs) < 2
                    eps = kwargs[1]
                    mu = 1.0
                    nk = epston(eps)
                else
                    eps, mu = kwargs
                    nk = epston(eps*mu)
                end
            end
        end
        return new(id, eps, mu, nk)
    end
end
# Material can be specified as Material(eps, mu), or just Material(eps), this will by default take mu = 1.0

"""
    Interface
This is struct object that defines interface datastruct
# Arguments
- `mat1::Material`: 1st material
- `mat2::Material`: 2nd material
# Examples
```julia-repl
julia> M1 = Material("nk", 2, 3)
Material("nk", -5.0 + 12.0im, 1.0, 2.0 + 3.0im)
julia> M2 = Material("epsmu", 2+3im, 1.0 + 0.0im)
Material("epsmu", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im)
julia> I = Interface(M1, M2)
Interface(Material("nk", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), Material("epsmu", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 0.0)
```
"""
mutable struct Interface
    mat1::Material
    mat2::Material
    sig::Number
    function Interface(mat1 = Material(), mat2 = Material(), args...)
        if isempty(args)
            sig = 0.0
        else
            sig = args[1]
        end
        return new(mat1, mat2, sig)
     end
end

"""
    Layer
This is struct object that defines layer datastruct
# Arguments
- `mat::Material`: Layer material
- `d::Real`: Layer thickness
# Examples
```julia-repl
julia> L = Layer(M1, 200)
Layer(Material("nk", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 200)
```
"""
mutable struct Layer
    mat::Material
    d::Real
    function Layer(mat = Material(), kwargs...)
        if isempty(kwargs)
            d = 0.0
        else
            d = kwargs[1]
        end
        return new(mat, d)
    end
end

function extract_params(input::Union{Material, Interface, Layer})
    if typeof(input) == Material
        return input.eps, input.mu
    elseif typeof(input) == Interface
        return input.mat1, input.mat2, input.sig
    else
        return input.id, input.mat, input.d
    end
end

function update_interface_list(L::Vector{Layer}, sig::Vector{<:Number}, nl::Integer)
    Is =  Vector{Interface}(undef, nl -1)
    for i = 1:(nl -1)
        Is[i] = Interface(L[i].mat, L[i+1].mat, sig[i])
    end
    return Is
end

"""
    Stack
This is struct object that defines layer datastruct
# Arguments
- `Layers::Vector{Layer}`: List of layers forming the stack. Incidence light is assumed to be from the leftmost layer.
- `Sigmas::Vector{<:Number}`: List of conductivities at the interfaces of the layers.
# Examples
```julia-repl
julia> L0 = Layer()
Layer(Material("nk", 1.0, 1.0, 1.0), 0.0)

julia> L1 = Layer(M1, 100)
Layer(Material("nk", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 100)

julia> L2 = Layer(M2, 600)
Layer(Material("epsmu", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 600)

julia> S = Stack([L0, L1, L2, L0], [0, 0, 0.2im])
Stack(Layer[Layer(Material("nk", 1.0, 1.0, 1.0), 0.0), Layer(Material("nk", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 100), Layer(Material("epsmu", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 600), Layer(Material("nk", 1.0, 1.0, 1.0), 0.0)], ComplexF64[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.2im], Interface[Interface(Material("nk", 1.0, 1.0, 1.0), Material("nk", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 0.0 + 0.0im), Interface(Material("nk", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), Material("epsmu", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 0.0 + 0.0im), Interface(Material("epsmu", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), Material("nk", 1.0, 1.0, 1.0), 0.0 + 0.2im)])
```
"""
mutable struct Stack
    Layers::Vector{Layer} #List of layers
    Sigmas::Vector{<:Number}
    Interfs::Vector{Interface}  # List of Interfaces
    function Stack(Layers, Sigmas)
        nL = length(Layers)
        Interfs = update_interface_list(Layers, Sigmas, nL)
        return new(Layers, Sigmas, Interfs)
    end
end

function modify_stack_layers(S::Stack, r::UnitRange{<:Integer}, L = Layer[], sigs = Number[])
    #This isn't finished yet, need to modified to add the interface update as well
    splice!(S.Layers, r, L)
    splice!(S.Sigmas, r, sigs)
    S.Interfs = update_interface_list(S.Layers, S.Sigmas, length(S.Layers))
end

"""
    revert_stack(S::Stack)
This is a function that reverts the order of the layers in the stuck
# Arguments
- `S::Stack`: Original stack
"""
function revert_stack(S::Stack)
    Snew = Stack(S.Layers[end:-1:1], S.Sigmas[end:-1:1])
    return Snew
end

function betz(mat::Material, k::Number, q::Real)
    return sqrt(mat.eps*mat.mu*k^2 - q^2 + 0.0im)
end

"""
    intface_rt(p::Integer, Iij::Interface, k0::Real, kp::Vector{<:Real})
This function calculates the Fresnel interface reflection and transmission coeficients
# Arguments
- `p::Integer`: Polarization of the incident plane wave, `p` = 1 corresponds to "s"-polarization, `p` = 2 corresponds to "p"-polarization.
- `Iij::Interface`: Interface object
- `k0::Real`: Free space light wavenumber
- `kp::Vector{<:Real}`: In-plane wavevector.
"""
function intface_rt(p::Integer, Iij::Interface, k0::Real, kp::Vector{<:Real})
    # Nottice: This is also coherent with Novotny-Hecht, an r, t
    #in both casesare just the rations of the electric fields, in contrast to formulations with magnetic field
    @assert p == 1 || p == 2
    local q = norm(kp)
    mat1, mat2, σ = extract_params(Iij)
    kz1 = betz(mat1, k0, q)
    kz2 = betz(mat2, k0, q)
    gs = mat1.mu/mat2.mu
    if p == 1 # s-polarization
        g = gs
    else
        g = mat1.eps/mat2.eps
    end
    r = (kz1 - g*kz2)/(kz1 + g*kz2)
    t = sqrt(g/gs)*(1 + r)
    return r, t
end

function Matrix2x2(r, t, d)
    local pp = cis(d)
    local cp = cis(-d)
    return [cp r*cp; r*pp pp]./t
end

function Matrix2x2_inc(r, t, d)
end

"""
    tmm_matrix(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
This function calculates the Fresnel reflection ``r_{\\rm s, p}`` and transmission ``t_{\\rm s, p}`` coeficients for layer stack using transfer matrix method.
# Arguments
- `p::Integer`: Polarization of the incident plane wave, `p` = 1 corresponds to "s"-polarization, `p` = 2 corresponds to "p"-polarization.
- `lambda::Real`: Free space light wavelength
- `kp::Vector{<:Real}`: In-plane wavevector.
- `S::Stack`: Layer stack fpr which the reflection and transmission coeficients will be caclulated. The light is assumed to be incident from the leftmost layer.
"""
function tmm_matrix(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
    local k0 = 2.0*pi/lambda;
    local q = sqrt(kp[1]^2 + kp[2]^2);
    Mglobal = [complex(1.0, 0.0) complex(0.0, 0.0); complex(0.0, 0.0) complex(1.0, 0.0)]
    Mlocal = zeros(ComplexF64, 2, 2)
    L::Layer = Layer()
    I::Interface = Interface()
    nL = length(S.Layers)
    for i = 1:(nL-1)
        L = S.Layers[i]
        I = S.Interfs[i]
        r, t = intface_rt(p, I, k0, kp)
        di = betz(L.mat, k0, q)*L.d
        Mlocal = Matrix2x2(r, t, di)
        Mglobal  = Mglobal*Mlocal
    end
    r = Mglobal[2,1]/Mglobal[1,1]
    t = 1.0/Mglobal[1,1]
    return r, t
end

"""
    RT_calc(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
This function calculates the reflectance `R` and transmittance `T` for layer stack using transfer matrix method.
# Arguments
- `p::Integer`: Polarization of the incident plane wave, `p` = 1 corresponds to `s`-polarization, `p` = 2 corresponds to `p`-polarization.
- `lambda::Real`: Free space light wavelength
- `kp::Vector{<:Real}`: In-plane wavevector.
- `S::Stack`: Layer stack fpr which the reflection and transmission coeficients will be caclulated.
# Example
Below we will use this code to calculate reflection coeficients of a layered strucutres consisting of 200 nm thick layer with constant dielectric function ``\\varepsilon_2 = 12.0 + 0.1{\\rm i}`` on an infinite substrate with refractive index ``n_3 = 1.45``. The superstrate is asumed to be air with refractive index ``n_1 = 1``.
## Example Cose
```julia
    deg = pi/180  # degree to radian conversion
    θ = 25*deg    # azimuthal angle of the incidence plane wave
    ϕ = 0*deg     # polar angle of the incidence plane wave

    mat1 = Material()  # Air material is initiated with "Material" object
                       # with empty argument which by default is vacuum.
    mat2 = Material("epsmu", 12 + 0.1im, 1.0)    # material 2 with ``\\varepsilon = 12 + 0.1im`` and ``\\mu = 1``.
    mat3 = Material("nk", 1.45, 0.0)             # material 3 is with constant refractive index n = 1.45, k = 0.
    L1 = Layer(mat1, 0)             # layer 1, is air with 0 thickness which the code will interprate as infinite medium.
    L2 = Layer(mat2, 200)           # layer 2 with thickness = 200 nm
    L3 = Layer(mat3, 0)             # layer 3 inifinite substrate
    S = Stack([L1, L2, L3], zeros(2))        # layer stack. The incidence light comes from the left.
    lambda = range(300, 900, length = 500)   # wavelength array
    Rs, Ts = zeros(500), zeros(500)          # allocated arrays for 's'-polarized R, T
    Rp, Tp = zeros(500), zeros(500)          # allocated arrays for 'p'-polarized R, T

    # Now we start the wavelength sweeping loop
    for (i, λ) = enumerate(lambda)
        k0 = 2*π/λ                 # wavenumber in free space
        kx = k0*sin(θ)*cos(ϕ)      # x component of the in-plane wavevector of the incident light
        ky = k0*sin(θ)*sin(ϕ)      # y component of the in-plane wavevector of the incident light
        # note that we assumed that the wave propagates along the 'z'-axis.
        Rs[i], Ts[i] = RT_calc(1, λ, [kx, ky], S)
        Rp[i], Tp[i] = RT_calc(2, λ, [kx, ky], S)
    end
    fig, ax = subplots(1,1)
    ax[:plot](lambda, Rs, label = "R\$_{s}\$", color = "black", linestyle = "-")
    ax[:plot](lambda, Ts, label = "T\$_{s}\$", color = "black", linestyle = "--")
    ax[:plot](lambda, Rp, label = "R\$_{p}\$", color = "red", linestyle = "-")
    ax[:plot](lambda, Tp, label = "T\$_{p}\$", color = "red", linestyle = "--")
    ax[:set_xlabel]("Wavelength (nm)")
    ax[:set_ylabel]("Reflection/Transmission")
    ax[:legend](frameon = false)
```⠀
"""
function RT_calc(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)
    local k0 = 2.0*pi/lambda;
    local q = sqrt(kp[1]^2 + kp[2]^2);
    r, t = tmm_matrix(p, lambda, kp, S)
    mat1 = S.Layers[1].mat
    matN  = S.Layers[end].mat
    kz1 = betz(mat1, k0, q)
    kzN = betz(matN, k0, q)
    R = abs2(r)
    T = abs2(t)
    T = (mat1.mu/matN.mu)*real(kzN/kz1)*T
    return R, T
end

"""
    ellipso_tmm(lambda::Real, kp::Vector{<:Real}, S::Stack)
This function returns the ellipsometric paramets of layer stack. These parameters are calculated as ``{\\rm tan}(\\psi)=\\vert r_p/r_s\\vert`` and
``{\\rm cos}(\\Delta) = {\\rm cos}(\\phi_p - \\phi_s)``, where ``\\phi_{s,p} = {\\rm angle}(r_{s,p})``
# Arguments
- `lambda::Real`: Free space light wavelength
- `kp::Vector{<:Real}`: In-plane wavevector.
- `S::Stack`: Layer stack fpr which the reflection and transmission coeficients will be caclulated.
"""
function ellipso_tmm(lambda::Real, kp::Vector{<:Real}, S::Stack)
    rs, ts =  tmm_matrix(1, lambda, kp, S)
    rp, tp =  tmm_matrix(2, lambda, kp, S)
    rr = rp/rs
    return abs(rr), cos(angle(rr))
end

#=======================  Example of usage ====================================
using PyPlot
pygui(true)

mat1 = Material()
mat2 = Material("epsmu", 12 + 0.1im)
mat3 = Material("nk", 1.45, 0.0)
L1 = Layer(mat1, 0)
L2 = Layer(mat2, 200)
L3 = Layer(mat3, 0)
S = Stack([L1, L2, L3], zeros(2))

deg = pi/180
θ = 25*deg
ϕ = 0*deg
lambda = range(300, 900, length = 500)
Rs, Ts = zeros(500), zeros(500)
Rp, Tp = zeros(500), zeros(500)
for (i, λ) = enumerate(lambda)
    k0 = 2*π/λ
    kx = k0*sin(θ)*cos(ϕ)
    ky = k0*sin(θ)*sin(ϕ)
    Rs[i], Ts[i] = RT_calc(1, λ, [kx, ky], S)
    Rp[i], Tp[i] = RT_calc(2, λ, [kx, ky], S)
end

fig, ax = subplots(1,1, figsize = (6, 4))
ax[:plot](lambda, Rs, label = "R\$_{s}\$", color = "black", linestyle = "-")
ax[:plot](lambda, Ts, label = "T\$_{s}\$", color = "black", linestyle = "--")
ax[:plot](lambda, Rp, label = "R\$_{p}\$", color = "red", linestyle = "-")
ax[:plot](lambda, Tp, label = "T\$_{p}\$", color = "red", linestyle = "--")
ax[:set_xlabel]("Wavelength (nm)")
ax[:set_ylabel]("Reflection/Transmission")
ax[:legend](frameon = false)
tight_layout()
savefig("C:/Users/vmkhi/Documents/Github/DipoleLattice/docs/Example Pictures/ExampleTMM.png", dpi = 600, transparent = true)

=#
