var documenterSearchIndex = {"docs":
[{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  Lattice2D","category":"page"},{"location":"Lattice/#DipoleLattice.Lattice2D","page":"Lattice","title":"DipoleLattice.Lattice2D","text":"lattice2D\n\nThis data struct constructs 2D lattice object and outputs its properties.\n\nArguments\n\na::Real: Length of lattice vector a1, this is assumed to be along x-axis.\nb::Real : Length of lattice vector a2\nψ::Real : Angle between lattice vector a1 and a2.\n\n\n\n\n\n","category":"type"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"(Image: Output)","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"<img src=\"Pictures/lattice2D.png\" width=\"75%\"/>","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  ConstructWZC(s1::Vector{<:Real}, s2::Vector{<:Real})","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  PlotBZ(L::Lattice2D)","category":"page"},{"location":"Lattice/#DipoleLattice.PlotBZ-Tuple{Lattice2D}","page":"Lattice","title":"DipoleLattice.PlotBZ","text":"PlotBZ(L::Lattice2D)\n\nArguments\n\nL::Lattice2D: 2D lattice object\n\n\n\n\n\n","category":"method"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"  rt_P_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αe::Array{<:Number, 2})","category":"page"},{"location":"#DipoleLattice.rt_P_dipole-Tuple{String, Real, Real, Real, Vector{<:Real}, Lattice2D, Matrix{<:Number}}","page":"Dipole Lattice","title":"DipoleLattice.rt_P_dipole","text":"rt_P_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αe::Array{<:Number, 2})\n\nThis function returns the field rflection and transmission coeficienets of a 2D array of dipoles with only electric dipole moment.\n\nArguments\n\npol::String: Polarization of incidence plane wave\nλ::Real : The free space wavelength of the incidence light\nθ::Real : The azimuthal angle of the incidence plane wave\nϕ::Real : The polar angle of the the incidence plane wave\nϵμh::Vector{<:Real} : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].\nlattice::Lattice2D : Lattice objects that defines teh properties of 2D-lattice\nαe : 3x3 electric polarizability tensor.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"rt_M_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αm::Array{<:Number, 2})","category":"page"},{"location":"#DipoleLattice.rt_M_dipole-Tuple{String, Real, Real, Real, Vector{<:Real}, Lattice2D, Matrix{<:Number}}","page":"Dipole Lattice","title":"DipoleLattice.rt_M_dipole","text":"rt_M_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αm::Array{<:Number, 2})\n\nThis function returns the field rflection and transmission coeficienets of a 2D array of dipoles with only magnetic dipole moments.\n\nArguments\n\npol::String: Polarization of incidence plane wave\nλ::Real : The free space wavelength of the incidence light\nθ::Real : The azimuthal angle of the incidence plane wave\nϕ::Real : The polar angle of the the incidence plane wave\nϵμh::Vector{<:Real} : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].\nlattice::Lattice2D : Lattice objects that defines teh properties of 2D-lattice\nαm : 3x3 magnetic polarizability tensor.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"  rt_PM_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αem::Array{<:Number, 2})","category":"page"},{"location":"#DipoleLattice.rt_PM_dipole-Tuple{String, Real, Real, Real, Vector{<:Real}, Lattice2D, Matrix{<:Number}}","page":"Dipole Lattice","title":"DipoleLattice.rt_PM_dipole","text":"rt_PM_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::lattice2D, αem::Array{<:Number, 2})\n\nThis function returns the field rflection and transmission coeficienets of a 2D array of dipoles with both electric and magnetic dipole moments.\n\nArguments\n\npol::String: Polarization of incidence plane wave\nλ::Real : The free space wavelength of the incidence light\nθ::Real : The azimuthal angle of the incidence plane wave\nϕ::Real : The polar angle of the the incidence plane wave\nϵμh::Vector{<:Real} : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].\nlattice::Lattice2D : Lattice objects that defines teh properties of 2D-lattice\nαem : 6x6 magneto-electrical polarizability tensor.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"  DyadSum(k::Number, rv::Vector{<:Real}, kp::Vector{<:Real}, lattice::Lattice2D, GorH::String; Nmax::Integer = 5)","category":"page"},{"location":"#DipoleLattice.DyadSum-Tuple{Number, Vector{<:Real}, Vector{<:Real}, Lattice2D, String}","page":"Dipole Lattice","title":"DipoleLattice.DyadSum","text":"DyadSum(k::Number, rv::Vector{<:Real}, kp::Vector{<:Real}, L::Lattice2D, GorH::String; Nmax::Integer = 5)\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  epston(epsilon::Number)","category":"page"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  Material","category":"page"},{"location":"tmm/#DipoleLattice.Material","page":"TMM","title":"DipoleLattice.Material","text":"Material\n\nThis is struct object that defines material datastruct\n\nArguments\n\nid::String: id=\"nk\" or id = \"epsmu\" defines the input type. If the id=\"nk\", we need to specify the refractive index. If id = \"epsmu\", the material will be defined through dielectric function and magnetic permeability.\neps::Number: Dielectric function of the material\nmu::Number: magnetic permeability of the material\nnk::Number: complex refractive index of the material\n\nExamples\n\njulia> M1 = Material(\"nk\", 2, 3)\nMaterial(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im)\njulia> M2 = Material(\"epsmu\", 2+3im, 1.0 + 0.0im)\nMaterial(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im)\n\n\n\n\n\n","category":"type"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  Interface","category":"page"},{"location":"tmm/#DipoleLattice.Interface","page":"TMM","title":"DipoleLattice.Interface","text":"Interface\n\nThis is struct object that defines interface datastruct\n\nArguments\n\nmat1::Material: 1st material\nmat2::Material: 2nd material\n\nExamples\n\njulia> M1 = Material(\"nk\", 2, 3)\nMaterial(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im)\njulia> M2 = Material(\"epsmu\", 2+3im, 1.0 + 0.0im)\nMaterial(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im)\njulia> I = Interface(M1, M2)\nInterface(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 0.0)\n\n\n\n\n\n","category":"type"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  Layer","category":"page"},{"location":"tmm/#DipoleLattice.Layer","page":"TMM","title":"DipoleLattice.Layer","text":"Layer\n\nThis is struct object that defines layer datastruct\n\nArguments\n\nmat::Material: Layer material\nd::Real: Layer thickness\n\nExamples\n\njulia> L = Layer(M1, 200)\nLayer(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 200)\n\n\n\n\n\n","category":"type"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  Stack","category":"page"},{"location":"tmm/#DipoleLattice.Stack","page":"TMM","title":"DipoleLattice.Stack","text":"Stack\n\nThis is struct object that defines layer datastruct\n\nArguments\n\nLayers::Vector{Layer}: List of layers forming the stack. Incidence light is assumed to be from the leftmost layer.\nSigmas::Vector{<:Number}: List of conductivities at the interfaces of the layers.\n\nExamples\n\njulia> L0 = Layer()\nLayer(Material(\"nk\", 1.0, 1.0, 1.0), 0.0)\n\njulia> L1 = Layer(M1, 100)\nLayer(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 100)\n\njulia> L2 = Layer(M2, 600)\nLayer(Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 600)\n\njulia> S = Stack([L0, L1, L2, L0], [0, 0, 0.2im])\nStack(Layer[Layer(Material(\"nk\", 1.0, 1.0, 1.0), 0.0), Layer(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 100), Layer(Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 600), Layer(Material(\"nk\", 1.0, 1.0, 1.0), 0.0)], ComplexF64[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.2im], Interface[Interface(Material(\"nk\", 1.0, 1.0, 1.0), Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 0.0 + 0.0im), Interface(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 0.0 + 0.0im), Interface(Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), Material(\"nk\", 1.0, 1.0, 1.0), 0.0 + 0.2im)])\n\n\n\n\n\n","category":"type"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  revert_stack(S::Stack)","category":"page"},{"location":"tmm/#DipoleLattice.revert_stack-Tuple{Stack}","page":"TMM","title":"DipoleLattice.revert_stack","text":"revert_stack(S::Stack)\n\nThis is a function that reverts the order of the layers in the stuck\n\nArguments\n\nS::Stack: Original stack\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  intface_rt(p::Integer, Iij::Interface, k0::Real, kp::Vector{<:Real})","category":"page"},{"location":"tmm/#DipoleLattice.intface_rt-Tuple{Integer, Interface, Real, Vector{<:Real}}","page":"TMM","title":"DipoleLattice.intface_rt","text":"intface_rt(p::Integer, Iij::Interface, k0::Real, kp::Vector{<:Real})\n\nThis function calculates the Fresnel interface reflection and transmission coeficients\n\nArguments\n\np::Integer: Polarization of the incident plane wave, p = 1 corresponds to \"s\"-polarization, p = 2 corresponds to \"p\"-polarization.\nIij::Interface: Interface object\nk0::Real: Free space light wavenumber\nkp::Vector{<:Real}: In-plane wavevector.\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  tmm_matrix(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)","category":"page"},{"location":"tmm/#DipoleLattice.tmm_matrix-Tuple{Integer, Real, Vector{<:Real}, Stack}","page":"TMM","title":"DipoleLattice.tmm_matrix","text":"tmm_matrix(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)\n\nThis function calculates the Fresnel reflection r_rm s p and transmission t_rm s p coeficients for layer stack using transfer matrix method.\n\nArguments\n\np::Integer: Polarization of the incident plane wave, p = 1 corresponds to \"s\"-polarization, p = 2 corresponds to \"p\"-polarization.\nlambda::Real: Free space light wavelength\nkp::Vector{<:Real}: In-plane wavevector.\nS::Stack: Layer stack fpr which the reflection and transmission coeficients will be caclulated. The light is assumed to be incident from the leftmost layer.\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  RT_calc(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)","category":"page"},{"location":"tmm/#DipoleLattice.RT_calc-Tuple{Integer, Real, Vector{<:Real}, Stack}","page":"TMM","title":"DipoleLattice.RT_calc","text":"RT_calc(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)\n\nThis function calculates the reflectance R and transmittance T for layer stack using transfer matrix method.\n\nArguments\n\np::Integer: Polarization of the incident plane wave, p = 1 corresponds to s-polarization, p = 2 corresponds to p-polarization.\nlambda::Real: Free space light wavelength\nkp::Vector{<:Real}: In-plane wavevector.\nS::Stack: Layer stack fpr which the reflection and transmission coeficients will be caclulated.\n\nExample\n\nBelow we will use this code to calculate reflection coeficients of a layered strucutres consisting of 200 nm thick layer with constant dielectric function varepsilon_2 = 120 + 01rm i on an infinite substrate with refractive index n_3 = 145. The superstrate is asumed to be air with refractive index n_1 = 1.\n\nExample Cose\n\n    deg = pi/180  # degree to radian conversion\n    θ = 25*deg    # azimuthal angle of the incidence plane wave\n    ϕ = 0*deg     # polar angle of the incidence plane wave\n\n    mat1 = Material()  # Air material is initiated with \"Material\" object\n                       # with empty argument which by default is vacuum.\n    mat2 = Material(\"epsmu\", 12 + 0.1im, 1.0)    # material 2 with ``\\varepsilon = 12 + 0.1im`` and ``\\mu = 1``.\n    mat3 = Material(\"nk\", 1.45, 0.0)             # material 3 is with constant refractive index n = 1.45, k = 0.\n    L1 = Layer(mat1, 0)             # layer 1, is air with 0 thickness which the code will interprate as infinite medium.\n    L2 = Layer(mat2, 200)           # layer 2 with thickness = 200 nm\n    L3 = Layer(mat3, 0)             # layer 3 inifinite substrate\n    S = Stack([L1, L2, L3], zeros(2))        # layer stack. The incidence light comes from the left.\n    lambda = range(300, 900, length = 500)   # wavelength array\n    Rs, Ts = zeros(500), zeros(500)          # allocated arrays for 's'-polarized R, T\n    Rp, Tp = zeros(500), zeros(500)          # allocated arrays for 'p'-polarized R, T\n\n    # Now we start the wavelength sweeping loop\n    for (i, λ) = enumerate(lambda)\n        k0 = 2*π/λ                 # wavenumber in free space\n        kx = k0*sin(θ)*cos(ϕ)      # x component of the in-plane wavevector of the incident light\n        ky = k0*sin(θ)*sin(ϕ)      # y component of the in-plane wavevector of the incident light\n        # note that we assumed that the wave propagates along the 'z'-axis.\n        Rs[i], Ts[i] = RT_calc(1, λ, [kx, ky], S)\n        Rp[i], Tp[i] = RT_calc(2, λ, [kx, ky], S)\n    end\n    fig, ax = subplots(1,1)\n    ax[:plot](lambda, Rs, label = \"R$_{s}$\", color = \"black\", linestyle = \"-\")\n    ax[:plot](lambda, Ts, label = \"T$_{s}$\", color = \"black\", linestyle = \"--\")\n    ax[:plot](lambda, Rp, label = \"R$_{p}$\", color = \"red\", linestyle = \"-\")\n    ax[:plot](lambda, Tp, label = \"T$_{p}$\", color = \"red\", linestyle = \"--\")\n    ax[:set_xlabel](\"Wavelength (nm)\")\n    ax[:set_ylabel](\"Reflection/Transmission\")\n    ax[:legend](frameon = false)\n\n⠀\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  <img src=\"Pictures/exampleTMM.png\" width=\"60%\"/>","category":"page"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  ellipso_tmm(lambda::Real, kp::Vector{<:Real}, S::Stack)","category":"page"},{"location":"tmm/#DipoleLattice.ellipso_tmm-Tuple{Real, Vector{<:Real}, Stack}","page":"TMM","title":"DipoleLattice.ellipso_tmm","text":"ellipso_tmm(lambda::Real, kp::Vector{<:Real}, S::Stack)\n\nThis function returns the ellipsometric paramets of layer stack. These parameters are calculated as rm tan(psi)=vert r_pr_svert and rm cos(Delta) = rm cos(phi_p - phi_s), where phi_sp = rm angle(r_sp)\n\nArguments\n\nlambda::Real: Free space light wavelength\nkp::Vector{<:Real}: In-plane wavevector.\nS::Stack: Layer stack fpr which the reflection and transmission coeficients will be caclulated.\n\n\n\n\n\n","category":"method"}]
}
