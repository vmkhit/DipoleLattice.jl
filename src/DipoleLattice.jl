"""
    Ddipole Lattice
    Author: Vahagn Mkhitaryan, vmkhit@gmail.com
    Date: 2022/03/22
    Notes: This package provides set of functions to calculate the scattering of
           2D array of dipole scattereres with full 6x6 dipole polarizabilities.
"""
module DipoleLattice
    import SpecialFunctions
    import LinearAlgebra
    import DelimitedFiles
    import Dierckx
    using PyPlot
    include("Lattice.jl")
    #include("Wijers.jl")
    include("DyadSum.jl")
    include("RDyadSum.jl")
    include("TMM.jl")
    include("MulticomponentLattice.jl")
    include("Lattice_reflection.jl")
    include("../utils/utils.jl")
    export Lattice2D, ConstructWZC, PlotBZ
    export DyadSum, ScalarGreenFn, DyadGreenFn, RDyadSum
    export Material, Interface, Layer, Stack, revert_stack, intface_rt, tmm_matrix, RT_calc, ellipso_tmm
    export CircleIntersect, SphereIntersect, CheckClusterIntersect, PlotSpheres, alpha_disc, eps_Drude
    export alpha_sp, al_em
    export rt_P_dipole, rt_M_dipole, rt_PM_dipole, rt_PM_dipole_FullMatrix
    export get_propagating_gvectors, Make_alpha3x3, Make_alpha6x6, Dipoles, Spheres, getpolarizationvectors
    export RT_PM_MultiSpheres_v1, RT_PM_MultiSpheres_v2, RT_PM_DipoleCluster, RT_P_DipoleCluster
    export RT_M_DipoleCluster, EHFAR_PM_DipoleCluster, RT_PM_DipoleCluster_diffractive
end # module
