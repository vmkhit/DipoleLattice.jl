include("../Mie.jl")
using PyPlot
pygui(true)
lam = range(100, 5000, length = 1000)
ra = 300
eps_SiO2 = 3.44^2


ale = Array{ComplexF64}(undef, 1000)
alm = Array{ComplexF64}(undef, 1000)
for (i, li) = enumerate(lam)
    ale[i], alm[i] = alpha_sp(li, [1. , 1.], [eps_SiO2, 1.0], ra)
end

begin
    plot(lam, real(ale), label = "\${\\rm Re}\\{\\alpha_{\\rm ee}\\}\$")
    plot(lam, imag(ale), label = "\${\\rm Im}\\{\\alpha_{\\rm ee}\\}\$")
    plot(lam, real(alm), label = "\${\\rm Re}\\{\\alpha_{\\rm mm}\\}\$")
    plot(lam, imag(alm), label = "\${\\rm Im}\\{\\alpha_{\\rm mm}\\}\$")
    legend()
end
