var documenterSearchIndex = {"docs":
[{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  Lattice2D","category":"page"},{"location":"Lattice/#DipoleLattice.Lattice2D","page":"Lattice","title":"DipoleLattice.Lattice2D","text":"This data struct constructs 2D lattice object and outputs its properties.\n\nArguments\n\nid::String: String that specify one of the five 2D bravis lattices or a gneral 2D lattice as described below.\n\n\n\n\n\n","category":"type"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"<img src=\"../Pictures/lattice2D.png\" width=\"75%\"/>","category":"page"},{"location":"Lattice/#Description","page":"Lattice","title":"Description","text":"","category":"section"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"The id argument takes the following values. The table below shows the id values and corresponding lattice types.","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"Lattice type id\nSquare \"S\"\nHexagonal \"H\"\nCentered Rectangular \"RC\"\nPrimitive Rectangular \"RP\"\nOblique \"O\"\nGeneral \"\"","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"If id = \"S\" or id = \"H\" one has to specify the lattice vector length as an optional argument. Similarly, for the case when id = \"RC\" or id = \"RP\" one has to specify the lengths of the lattice vectors along each orthogonal direction. For Oblique lattice, besides the lattice vector lengths one has to specify also the angle between these vectors as the 3rd optional argument. If id is specified as an empty string id = \"\" one has to specify the direct lattice vectors veca_1 and veca_2 as an optional arguments.","category":"page"},{"location":"Lattice/#Examples","page":"Lattice","title":"Examples","text":"","category":"section"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"To define a square lattice with unit vector length = a = 5 we can call.","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  julia> L = Lattice(\"S\", 5)","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"Hexagonal lattice with unit vector length  a = 4 can be defined as","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  julia> L = Lattice(\"S\", 4)","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"Rectangular lattices with unit vector lengths a = 2, b = 4 can be defined as.","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  julia> L = Lattice(\"RC\", a, b)","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"The primitive rectangular lattice can be constructed using:","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  julia> L = Lattice(\"RP\", a, b)","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"The oblique lattice takes additional argument psi the angle between the vectors, assuming veca_1 is along the x-axis.","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  julia> L = Lattice(\"O\", a, b, ψ)","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"Finally the general lattice with unit vectors veca_1 = a_1x a_1y and veca_1 = a_2x a_2y.","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  julia> L = Lattice(\"\", a1, a2)","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  ConstructWZC(R::Array{<:Real})","category":"page"},{"location":"Lattice/#DipoleLattice.ConstructWZC-Tuple{Array{<:Real}}","page":"Lattice","title":"DipoleLattice.ConstructWZC","text":"Arguments\n\nR::Array{<:Real}: Unit cell vectors of the 2D lattice, a_1 = R[1, :], a_2 = R[2, :].\n\n\n\n\n\nArguments\n\nR::Array{<:Real}: Unit cell vectors of the 2D lattice, a_1 = R[1, :], a_2 = R[2, :].\n\n\n\n\n\n","category":"method"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"This function implements the construction of the BZ of a 2D lattice following the procedure provided in reference Thompson, I., and Linton, C. M. (2010). \"Guided surface waves on one-and two-dimensional arrays of spheres\"., SIAM Journal on Applied Mathematics, 70(8), 2975-2995.","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  MakeKpath(verts::Array{<:Real, 2}, res::Union{Integer, Vector{<:Integer}}; close::Bool=true)","category":"page"},{"location":"Lattice/#DipoleLattice.MakeKpath-Tuple{Matrix{<:Real}, Union{Integer, Vector{<:Integer}}}","page":"Lattice","title":"DipoleLattice.MakeKpath","text":"Generates linear interpolation paths between the vertices provided in the list vertex.\n\nArguments\n\nverts::Array{<:Real, 2}:  Is a 2D arrar of points v_ix = verts[i, 1] and v_iy = verts[i, 2].\nres::Union{Integer, Vector{<:Integer}}: Number of sampling points between vertices. If it is specified as an integer N,   it will take a N sampling point between any pair of vertives. One can specify also list of integers with number   of points between each consequitive vertices.\nclose::Bool=true: Boolian that defines if the path is closed or not.\n\n\n\n\n\nGenerates linear interpolation paths between the vertices provided in the list vertex.\n\nArguments\n\nverts::Array{<:Real, 2}:  Is a 2D arrar of points v_ix = verts[i, 1] and v_iy = verts[i, 2].\nres::Union{Integer, Vector{<:Integer}}: Number of sampling points between vertices. If it is specified as an integer N,   it will take a N sampling point between any pair of vertives. One can specify also list of integers with number   of points between each consequitive vertices.\nclose::Bool=true: Boolian that defines if the path is closed or not.\n\n\n\n\n\n","category":"method"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  PointsInPolygon(v::Array{<:Real, 2}, n::Integer)","category":"page"},{"location":"Lattice/#DipoleLattice.PointsInPolygon-Tuple{Matrix{<:Real}, Integer}","page":"Lattice","title":"DipoleLattice.PointsInPolygon","text":"Generate points in a polygon. The algorithm devide the polygon into triangles by connecting the centroid of the polygon to each vertex.   Then each triangle is subdevided into smaller triangles by taking n points along each edge.\n\nArguments\n\nv::Array{<:Real, 2}: Vertices of the polygon, specified as v_ix = verts[i, 1] and v_iy = verts[i, 2].\nn::Integer: Number of subdivision along each edge of sub-triangles.\n\n\n\n\n\n","category":"method"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  PlotBZ(L::Lattice2D)","category":"page"},{"location":"Lattice/#DipoleLattice.PlotBZ-Tuple{Lattice2D}","page":"Lattice","title":"DipoleLattice.PlotBZ","text":"Arguments\n\nL::Lattice2D: 2D lattice object\n\nThis function plots the 1st Brillouin zone of the 2D lattice\n\n\n\n\n\n","category":"method"},{"location":"Lattice/#Example-BZ-construction","page":"Lattice","title":"Example BZ construction","text":"","category":"section"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"  julia> L = Lattice(\"O\", 1, 1, π/4), PlotBZ(L)","category":"page"},{"location":"Lattice/","page":"Lattice","title":"Lattice","text":"<img src=\"../Pictures/PlotBZ.png\" width=\"75%\"/>","category":"page"},{"location":"rt/","page":"Reflection and Transmission","title":"Reflection and Transmission","text":"  rt_P_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αe::Array{<:Number, 2})","category":"page"},{"location":"rt/#DipoleLattice.rt_P_dipole-Tuple{String, Real, Real, Real, Vector{<:Real}, Lattice2D, Matrix{<:Number}}","page":"Reflection and Transmission","title":"DipoleLattice.rt_P_dipole","text":"rt_P_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αe::Array{<:Number, 2})\n\nThis function returns the field rflection and transmission coeficienets of a 2D array of dipoles with only electric dipole moment.\n\nArguments\n\npol::String: Polarization of incidence plane wave\nλ::Real : The free space wavelength of the incidence light\nθ::Real : The azimuthal angle of the incidence plane wave\nϕ::Real : The polar angle of the the incidence plane wave\nϵμh::Vector{<:Real} : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].\nlattice::Lattice2D : Lattice objects that defines teh properties of 2D-lattice\nαe : 3x3 electric polarizability tensor.\n\n\n\n\n\n","category":"method"},{"location":"rt/","page":"Reflection and Transmission","title":"Reflection and Transmission","text":"rt_M_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αm::Array{<:Number, 2})","category":"page"},{"location":"rt/#DipoleLattice.rt_M_dipole-Tuple{String, Real, Real, Real, Vector{<:Real}, Lattice2D, Matrix{<:Number}}","page":"Reflection and Transmission","title":"DipoleLattice.rt_M_dipole","text":"rt_M_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αm::Array{<:Number, 2})\n\nThis function returns the field rflection and transmission coeficienets of a 2D array of dipoles with only magnetic dipole moments.\n\nArguments\n\npol::String: Polarization of incidence plane wave\nλ::Real : The free space wavelength of the incidence light\nθ::Real : The azimuthal angle of the incidence plane wave\nϕ::Real : The polar angle of the the incidence plane wave\nϵμh::Vector{<:Real} : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].\nlattice::Lattice2D : Lattice objects that defines teh properties of 2D-lattice\nαm : 3x3 magnetic polarizability tensor.\n\n\n\n\n\n","category":"method"},{"location":"rt/","page":"Reflection and Transmission","title":"Reflection and Transmission","text":"  rt_PM_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::Lattice2D, αem::Array{<:Number, 2})","category":"page"},{"location":"rt/#DipoleLattice.rt_PM_dipole-Tuple{String, Real, Real, Real, Vector{<:Real}, Lattice2D, Matrix{<:Number}}","page":"Reflection and Transmission","title":"DipoleLattice.rt_PM_dipole","text":"rt_PM_dipole(pol::String, λ::Real, θ::Real, ϕ::Real, ɛμh::Vector{<:Real}, lattice::lattice2D, αem::Array{<:Number, 2})\n\nThis function returns the field rflection and transmission coeficienets of a 2D array of dipoles with both electric and magnetic dipole moments.\n\nArguments\n\npol::String: Polarization of incidence plane wave\nλ::Real : The free space wavelength of the incidence light\nθ::Real : The azimuthal angle of the incidence plane wave\nϕ::Real : The polar angle of the the incidence plane wave\nϵμh::Vector{<:Real} : Dielectric function and magnetic permitivities of the host medium. We use  ϵμh = [ϵh, μh].\nlattice::Lattice2D : Lattice objects that defines teh properties of 2D-lattice\nαem : 6x6 magneto-electrical polarizability tensor.\n\n\n\n\n\n","category":"method"},{"location":"#DipoleLattice.jl","page":"Dipole Lattice","title":"DipoleLattice.jl","text":"","category":"section"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"This package provides set of routines to calculate scattering properties of free standing and supported 2D periodic arrays of magnetoelectric dipolar particles. It is assumed that the scattering properties of these particles can be given by 6x6 dipole polarizability tensor alpha","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"alpha = beginbmatrix\nalpha_ee  alpha_em\nalpha_me  alpha_mm\nendbmatrix","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"where, each of alpha_i j, i j = e m are three by three matrices that describe the particle electric and magnetic polarizabilities as well as their cross-coupling, such that the for a given external fields mathbfE^rm ext and mathbfH^rm ext the induced electric and magnetic dipole moments are given by","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"beginbmatrix\nmathbfp\nmathbfm\nendbmatrix = beginbmatrix\nalpha_ee  alpha_em\nalpha_me  alpha_mm\nendbmatrix\nbeginbmatrix\nmathbfE^rm ext \nmathbfH^rm ext\nendbmatrix","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"The electric and magnetic fields at position mathbfr produced by such an magnetoelectric scatterer located at some position mathbfr_0 can be given using the so called electric and magnetic dyadic green's functions of the electromagnetic field","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"  beginbmatrix\n  mathbfE \n  mathbfH^prime\n  endbmatrix =\n  beginbmatrix\n    mathcalG^neq_E(mathbfr-mathbfr_0)  mathrmi k mathcalG^neq_H(mathbfr-mathbfr_0)\n    -mathrmi k mathcalG^neq_H(mathbfr-mathbfr_0)  mathcalG^neq_E(mathbfr-mathbfr_0)\n  endbmatrix\n  beginbmatrix\n  mathbfp^prime \n  mathbfm^prime\n  endbmatrix = mathcalG_rm 6x6(mathbfr-mathbfr_0)beginbmatrix\n  mathbfp^prime \n  mathbfm^prime\n  endbmatrix","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"where for the sake of symmetry the following normalized quantities are introduced mathbfp^prime = mathbfpvarepsilon, mathbfm^prime = Zmathbfm and mathbfH^prime = ZmathbfH, with varepsilon and mu being the dielectric function and permeability of the host medium and Z = sqrtmuvarepsilon. The dyadic Green functons are given as","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"  mathcalG_E^neq(mathbfr-mathbfr) = leftI k^2 +nablaotimesnablarightfracmathrme^mathrmi k vert mathbfr-mathbfr vertvertmathbfr-mathbfrvert\n  mathcalG_H^neq(mathbfr-mathbfr) = dfrac1k^2 nablatimes mathcalG^neq_E(mathbfr-mathbfr) = I  times nablafracmathrme^mathrmi kvertmathbfr-mathbfrvertvertmathbfr-mathbfrvert","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"with k=k_0sqrtvarepsilonmu, k_0 = omegac. For short notations we will introduce column vectors boldsymbolmu = mathbfp mathbfm^prime^T and mathbfF = mathbfE mathbfH^prime^T, where superscript T stands for transpose. With this notations, we can write a self-consistent system of equations for the dipoles in the array","category":"page"},{"location":"","page":"Dipole Lattice","title":"Dipole Lattice","text":"  boldsymbolmu_i^alpha = mathbfF^rm ext(mathbfr_ialpha) + sum_j neq i mathcalG_rm 6x6(mathbfr_ialpha-mathbfr_j alpha)boldsymbolmu_j^alpha + sum_jsum_beta neq alpha mathcalG_rm 6x6(mathbfr_ialpha-mathbfr_j beta)boldsymbolmu_j^beta","category":"page"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  epston(epsilon::Number)","category":"page"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  Material","category":"page"},{"location":"tmm/#DipoleLattice.Material","page":"TMM","title":"DipoleLattice.Material","text":"Material\n\nThis is struct object that defines material datastruct\n\nArguments\n\nid::String: id=\"nk\" or id = \"epsmu\" defines the input type. If the id=\"nk\", we need to specify the refractive index. If id = \"epsmu\", the material will be defined through dielectric function and magnetic permeability.\neps::Number: Dielectric function of the material\nmu::Number: magnetic permeability of the material\nnk::Number: complex refractive index of the material\n\nExamples\n\njulia> M1 = Material(\"nk\", 2, 3)\nMaterial(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im)\njulia> M2 = Material(\"epsmu\", 2+3im, 1.0 + 0.0im)\nMaterial(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im)\n\n\n\n\n\n","category":"type"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  Interface","category":"page"},{"location":"tmm/#DipoleLattice.Interface","page":"TMM","title":"DipoleLattice.Interface","text":"Interface\n\nThis is struct object that defines interface datastruct\n\nArguments\n\nmat1::Material: 1st material\nmat2::Material: 2nd material\n\nExamples\n\njulia> M1 = Material(\"nk\", 2, 3)\nMaterial(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im)\njulia> M2 = Material(\"epsmu\", 2+3im, 1.0 + 0.0im)\nMaterial(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im)\njulia> I = Interface(M1, M2)\nInterface(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 0.0)\n\n\n\n\n\n","category":"type"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  Layer","category":"page"},{"location":"tmm/#DipoleLattice.Layer","page":"TMM","title":"DipoleLattice.Layer","text":"Layer\n\nThis is struct object that defines layer datastruct\n\nArguments\n\nmat::Material: Layer material\nd::Real: Layer thickness\n\nExamples\n\njulia> L = Layer(M1, 200)\nLayer(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 200)\n\n\n\n\n\n","category":"type"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  Stack","category":"page"},{"location":"tmm/#DipoleLattice.Stack","page":"TMM","title":"DipoleLattice.Stack","text":"Stack\n\nThis is struct object that defines layer datastruct\n\nArguments\n\nLayers::Vector{Layer}: List of layers forming the stack. Incidence light is assumed to be from the leftmost layer.\nSigmas::Vector{<:Number}: List of conductivities at the interfaces of the layers.\n\nExamples\n\njulia> L0 = Layer()\nLayer(Material(\"nk\", 1.0, 1.0, 1.0), 0.0)\n\njulia> L1 = Layer(M1, 100)\nLayer(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 100)\n\njulia> L2 = Layer(M2, 600)\nLayer(Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 600)\n\njulia> S = Stack([L0, L1, L2, L0], [0, 0, 0.2im])\nStack(Layer[Layer(Material(\"nk\", 1.0, 1.0, 1.0), 0.0), Layer(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 100), Layer(Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 600), Layer(Material(\"nk\", 1.0, 1.0, 1.0), 0.0)], ComplexF64[0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.2im], Interface[Interface(Material(\"nk\", 1.0, 1.0, 1.0), Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), 0.0 + 0.0im), Interface(Material(\"nk\", -5.0 + 12.0im, 1.0, 2.0 + 3.0im), Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), 0.0 + 0.0im), Interface(Material(\"epsmu\", 2 + 3im, 1.0 + 0.0im, 1.6741492280355401 + 0.895977476129838im), Material(\"nk\", 1.0, 1.0, 1.0), 0.0 + 0.2im)])\n\n\n\n\n\n","category":"type"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  revert_stack(S::Stack)","category":"page"},{"location":"tmm/#DipoleLattice.revert_stack-Tuple{Stack}","page":"TMM","title":"DipoleLattice.revert_stack","text":"revert_stack(S::Stack)\n\nThis is a function that reverts the order of the layers in the stuck\n\nArguments\n\nS::Stack: Original stack\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  intface_rt(p::Integer, Iij::Interface, k0::Real, kp::Vector{<:Real})","category":"page"},{"location":"tmm/#DipoleLattice.intface_rt-Tuple{Integer, Interface, Real, Vector{<:Real}}","page":"TMM","title":"DipoleLattice.intface_rt","text":"intface_rt(p::Integer, Iij::Interface, k0::Real, kp::Vector{<:Real})\n\nThis function calculates the Fresnel interface reflection and transmission coeficients\n\nArguments\n\np::Integer: Polarization of the incident plane wave, p = 1 corresponds to \"s\"-polarization, p = 2 corresponds to \"p\"-polarization.\nIij::Interface: Interface object\nk0::Real: Free space light wavenumber\nkp::Vector{<:Real}: In-plane wavevector.\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  tmm_matrix(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)","category":"page"},{"location":"tmm/#DipoleLattice.tmm_matrix-Tuple{Integer, Real, Vector{<:Real}, Stack}","page":"TMM","title":"DipoleLattice.tmm_matrix","text":"tmm_matrix(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)\n\nThis function calculates the Fresnel reflection r_rm s p and transmission t_rm s p coeficients for layer stack using transfer matrix method.\n\nArguments\n\np::Integer: Polarization of the incident plane wave, p = 1 corresponds to \"s\"-polarization, p = 2 corresponds to \"p\"-polarization.\nlambda::Real: Free space light wavelength\nkp::Vector{<:Real}: In-plane wavevector.\nS::Stack: Layer stack fpr which the reflection and transmission coeficients will be caclulated. The light is assumed to be incident from the leftmost layer.\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  RT_calc(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)","category":"page"},{"location":"tmm/#DipoleLattice.RT_calc-Tuple{Integer, Real, Vector{<:Real}, Stack}","page":"TMM","title":"DipoleLattice.RT_calc","text":"RT_calc(p::Integer, lambda::Real, kp::Vector{<:Real}, S::Stack)\n\nThis function calculates the reflectance R and transmittance T for layer stack using transfer matrix method.\n\nArguments\n\np::Integer: Polarization of the incident plane wave, p = 1 corresponds to s-polarization, p = 2 corresponds to p-polarization.\nlambda::Real: Free space light wavelength\nkp::Vector{<:Real}: In-plane wavevector.\nS::Stack: Layer stack fpr which the reflection and transmission coeficients will be caclulated.\n\nExample\n\nBelow we will use this code to calculate reflection coeficients of a layered strucutres consisting of 200 nm thick layer with constant dielectric function varepsilon_2 = 120 + 01rm i on an infinite substrate with refractive index n_3 = 145. The superstrate is asumed to be air with refractive index n_1 = 1.\n\nExample Cose\n\n    deg = pi/180  # degree to radian conversion\n    θ = 25*deg    # azimuthal angle of the incidence plane wave\n    ϕ = 0*deg     # polar angle of the incidence plane wave\n\n    mat1 = Material()  # Air material is initiated with \"Material\" object\n                       # with empty argument which by default is vacuum.\n    mat2 = Material(\"epsmu\", 12 + 0.1im, 1.0)    # material 2 with ``\\varepsilon = 12 + 0.1im`` and ``\\mu = 1``.\n    mat3 = Material(\"nk\", 1.45, 0.0)             # material 3 is with constant refractive index n = 1.45, k = 0.\n    L1 = Layer(mat1, 0)             # layer 1, is air with 0 thickness which the code will interprate as infinite medium.\n    L2 = Layer(mat2, 200)           # layer 2 with thickness = 200 nm\n    L3 = Layer(mat3, 0)             # layer 3 inifinite substrate\n    S = Stack([L1, L2, L3], zeros(2))        # layer stack. The incidence light comes from the left.\n    lambda = range(300, 900, length = 500)   # wavelength array\n    Rs, Ts = zeros(500), zeros(500)          # allocated arrays for 's'-polarized R, T\n    Rp, Tp = zeros(500), zeros(500)          # allocated arrays for 'p'-polarized R, T\n\n    # Now we start the wavelength sweeping loop\n    for (i, λ) = enumerate(lambda)\n        k0 = 2*π/λ                 # wavenumber in free space\n        kx = k0*sin(θ)*cos(ϕ)      # x component of the in-plane wavevector of the incident light\n        ky = k0*sin(θ)*sin(ϕ)      # y component of the in-plane wavevector of the incident light\n        # note that we assumed that the wave propagates along the 'z'-axis.\n        Rs[i], Ts[i] = RT_calc(1, λ, [kx, ky], S)\n        Rp[i], Tp[i] = RT_calc(2, λ, [kx, ky], S)\n    end\n    fig, ax = subplots(1,1)\n    ax[:plot](lambda, Rs, label = \"R$_{s}$\", color = \"black\", linestyle = \"-\")\n    ax[:plot](lambda, Ts, label = \"T$_{s}$\", color = \"black\", linestyle = \"--\")\n    ax[:plot](lambda, Rp, label = \"R$_{p}$\", color = \"red\", linestyle = \"-\")\n    ax[:plot](lambda, Tp, label = \"T$_{p}$\", color = \"red\", linestyle = \"--\")\n    ax[:set_xlabel](\"Wavelength (nm)\")\n    ax[:set_ylabel](\"Reflection/Transmission\")\n    ax[:legend](frameon = false)\n\n⠀\n\n\n\n\n\n","category":"method"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  <img src=\"../Pictures/exampleTMM.png\" width=\"60%\"/>","category":"page"},{"location":"tmm/","page":"TMM","title":"TMM","text":"  ellipso_tmm(lambda::Real, kp::Vector{<:Real}, S::Stack)","category":"page"},{"location":"tmm/#DipoleLattice.ellipso_tmm-Tuple{Real, Vector{<:Real}, Stack}","page":"TMM","title":"DipoleLattice.ellipso_tmm","text":"ellipso_tmm(lambda::Real, kp::Vector{<:Real}, S::Stack)\n\nThis function returns the ellipsometric paramets of layer stack. These parameters are calculated as rm tan(psi)=vert r_pr_svert and rm cos(Delta) = rm cos(phi_p - phi_s), where phi_sp = rm angle(r_sp)\n\nArguments\n\nlambda::Real: Free space light wavelength\nkp::Vector{<:Real}: In-plane wavevector.\nS::Stack: Layer stack fpr which the reflection and transmission coeficients will be caclulated.\n\n\n\n\n\n","category":"method"}]
}
