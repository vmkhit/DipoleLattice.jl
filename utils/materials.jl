using SpecialFunctions
using LinearAlgebra
using DelimitedFiles
using Interpolations
using Dierckx
using DataFrames

eps_path = "C:\\Users\\vmkhitaryan\\Documents\\Github\\jcalc\\eps_data"

function importmat(filename, x, epsornk="n"; path = eps_path)
     file = filter(f-> startswith(f, filename) && endswith(f, ".txt"), readdir(eps_path))
     if !isempty(file)
         data = readdlm(joinpath(eps_path, file[1]))
         if size(data, 2) < 2
             throw("The data file must have at least two columns")
         elseif size(data, 2) == 2
             n_re = Spline1D(data[:,1], data[:, 2], k=1, bc="extrapolate")
             if epsornk=="n"
                 return complex.(n_re(x), 0.0)
             elseif epsornk == "e"
                 return complex.(n_re(x), 0.0).^2
             else
                 throw("The data type has to be 'e'-for eps, or 'n' for nk")
             end
         else
             n_re = Spline1D(data[:,1], data[:, 2], k=1, bc="extrapolate")
             n_im = Spline1D(data[:, 1], data[:, 3], k=1, bc="extrapolate")
             if epsornk == "n"
                 return complex.(n_re(x), n_im(x))
             elseif epsornk == "e"
                 return complex.(n_re(x), n_im(x)).^2
             else
                 throw("The data type has to be 'e'-for eps, or 'n' for nk")
             end
         end
     else
         throw("The specified file is not in the directory")
     end
end

"""
#=
Function Uasage
x = range(400, length = 100, stop = 5000)
@time nk = importmat("TiO2", x, "n")
@time ee = importmat("Ag_JC", x, "e")
plot(x, real(ee),  "blue", x, imag(ee), "black")
=#
"""

function eps2n(ee::ComplexF64)
    return sqrt(ee)
end
eps2n(ee::Real) = eps2n(complex(ee))

function Drude(w, eb, wp, gp)
    return  eb .- wp*wp./(w.*w .+ 1.0im.*w.*gp)
end

function eps_Au(w)
    hartree = 27.2116;              #  2 * Rydberg in eV
    tunit = 0.66 / hartree         #  time unit in fs
    rs = 3                   #  electron gas parameter
    eps0 = 10                #  background dielectric constant
    gammad = tunit/10         #  Drude relaxation rate
    #  density in atomic units
    density = 3/(4*pi *rs^3)
    #  plasmon energy
    wp = sqrt(4*pi*density)

    gammad = gammad*hartree
    wp = wp*hartree
    eps = eps0 - wp^2 ./(w.*(w + 1im*gammad))
    return wp, gammad,  eps0, eps
end

#wp, gam, epso, eps = eps_Au(0.5)

#wp, gam, epso, eps = eps_Au(0.5)

function eps_Drude_JC(w, name)
    names = ["Au", "Ag", "Cu"]
    wp   = [9.0833, 9.149, 8.82]
    gam0 = [0.071, 0.021, 0.0954]
    epsb = [9.1, 3.7, 4.5]
    ii =  findfirst(names, name)
    return epsb[ii] - wp[ii].^2 ./(w.^2 + 1im*w.*gam0[ii])
end


function create_df(name, models)
    table = DataFrame( Omegap = Number[],  gammap = Number[])
    push!(table, [0.5  0.01])
    test = "eps_data"
    if isdir("../$test") == 0
        mkdir("../$test")  # .. to move one folder up
    end
    writetable("Tests/$test/$name.csv",table)
end

function read_df(name)
    test = "eps_data"
    if isdir("../$test") == 0
        mkdir("../$test")  # .. to move one folder up
    end
    table = readtable("Tests/$test/$name.csv")
    return table
end



function epston(epsilon)
    n = sqrt(abs(epsilon) + real(epsilon))/sqrt(2)
    k = sqrt(abs(epsilon) - real(epsilon))/sqrt(2)
    return n+im*k
end

# Tauc-Lorentz Model: The imaginary part of this model constructed by multiplying the imaginary parts of the
# Tauc model and Lorentz model


function TL_N(ϵ_inf, E, E0_list, Eg, A_list, C_list)
  function TL(E, E0, Eg, A, C)
      if E > Eg
          ϵi = A*E0*C*(E-Eg)^2/(E*((E*E-E0*E0)^2+C*C*E*E))
      else
          ϵi = 0.0
      end

      aln = (Eg^2-E0^2)*E^2+Eg*Eg*C^2-(E0^2+3*Eg^2)*E0^2
      atg = (E^2-E0^2)*(E0^2+Eg^2)+Eg*Eg*C^2
      α = sqrt(complex(4.0*E0^2-C^2))
      γ = sqrt(complex(E0^2-0.5*C^2))
      ζ4 = (E^2-γ^2)^2+0.25*C*C*α^2

      ϵr =  (A*C*aln/(2*pi*ζ4*α*E0))*log((E0^2+Eg^2+α*Eg)/(E0^2+Eg^2-α*Eg))
      ϵr += (-A*atg/(pi*E0*ζ4))*(pi-atan((2*Eg+α)/C)+atan((α-2*Eg)/C))
      ϵr += (4*A*E0*Eg*(E^2-γ^2)/(pi*α*ζ4))*(atan((α+2*Eg)/C)+atan((α-2*Eg)/C))
      ϵr += (-A*E0*C*(E^2+Eg^2)/(pi*ζ4*E))*log(abs(E-Eg)/(E+Eg))
      ϵr += (2*A*E0*C*Eg/(pi*ζ4))*log(abs(E-Eg)*(E+Eg)/sqrt((E0^2-Eg^2)^2+C*C*Eg^2))

      return ϵr + 1.0im*ϵi
    end

    eps_TL =  ϵ_inf + 0.0im
    for i = 1:length(E0_list)
        eps_TL += TL(E, E0_list[i], Eg, A_list[i], C_list[i])
    end
    return eps_TL
end



function LDBB(metal, model, lambda0)
  c_const = 299792458.0
  ħ_const = 1.0545718001391127e-34
  e_const = 1.6021766208e-19

  # supported metals
  metals = ["Ag", "Au", "Cu", "Al", "Be", "Cr", "Ni", "Pd", "Pt", "Ti", "W"]
  lambda0 = lambda0*1e-9
  # supported models:
  # D - Drude
  # LD - Lorentz-Drude
  # BB - Brebhal
  models = ["D", "LD", "BB"]

  # material data for the LD model, from rakic et al
  LDdata = [
    [9.010 9.030 10.83 14.98 18.51 10.75 15.92 9.72 9.59 7.29 13.22];
    [0.845 0.760 0.575 0.523 0.084 0.168 0.096 0.330 0.333 0.148 0.206];
    [0.048 0.053 0.030 0.047 0.035 0.047 0.048 0.008 0.080 0.082 0.064];
    [0.065 0.024 0.061 0.227 0.031 0.151 0.100 0.649 0.191 0.899 0.054];
    [3.886 0.241 0.378 0.333 1.664 3.175 4.511 2.950 0.517 2.276 0.530];
    [0.816 0.415 0.291 0.162 0.100 0.121 0.174 0.336 0.780 0.777 1.004];
    [0.124 0.010 0.104 0.050 0.140 0.150 0.135 0.121 0.659 0.393 0.166];
    [0.452 0.345 1.056 0.312 3.395 1.305 1.334 0.555 1.838 2.518 1.281];
    [4.481 0.830 2.957 1.544 1.032 0.543 0.582 0.501 1.314 1.545 1.917];
    [0.011 0.071 0.723 0.166 0.530 1.149 0.106 0.638 0.547 0.187 0.706];
    [0.065 0.870 3.213 1.351 4.454 2.676 2.178 4.621 3.668 1.663 3.332];
    [8.185 2.969 5.300 1.808 3.183 1.970 1.597 1.659 3.141 2.509 3.580];
    [0.840 0.601 0.638 0.030 0.130 0.825 0.729 0.453 3.576 0.001 2.590];
    [0.916 2.494 4.305 3.382 1.802 1.335 6.292 3.236 8.517 1.762 5.836];
    [9.083 4.304 11.18 3.473 4.604 8.775 6.089 5.715 9.249 19.43 7.498];
    [5.646 4.384 0 0 0 0 0 0 0 0 0];
    [2.419 2.214 0 0 0 0 0 0 0 0 0];
    [20.29 13.32 0 0 0 0 0 0 0 0 0]
    ]

    # experimental data for the BB model, from rakic et al
  BBdata =[
    [9.01 9.03 10.83 14.98 18.51 10.75 15.92 9.72 9.59 7.29 13.22];
    [0.821 0.770 0.562 0.526 0.081 0.154 0.083 0.330 0.333 0.126 0.197];
    [0.049 0.050 0.030 0.047 0.035 0.048 0.022 0.009 0.080 0.067 0.057];
    [0.050 0.054 0.076 0.213 0.066 0.338 0.357 0.769 0.186 0.427 0.006];
    [0.189 0.074 0.056 0.312 2.956 4.256 2.820 2.343 0.498 1.877 3.689];
    [2.025 0.218 0.416 0.163 0.131 0.281 0.317 0.066 0.782 1.459 0.481];
    [1.894 0.742 0.562 0.013 0.277 0.115 0.606 0.694 0.031 0.463 3.754];
    [0.133 0.050 0.081 0.060 0.067 0.261 0.039 0.093 0.665 0.218 0.022];
    [0.067 0.035 0.047 0.315 3.962 3.957 0.120 0.497 1.851 0.100 0.277];
    [5.185 2.885 2.849 1.561 0.469 0.584 1.059 0.502 1.317 2.661 0.985];
    [0.665 0.349 0.469 0.042 3.167 0.252 1.454 0.027 0.096 0.506 0.059];
    [0.051 0.312 0.324 0.182 0.346 0.817 0.127 0.309 0.551 0.513 0.136];
    [0.019 0.083 0.113 1.587 2.398 2.218 1.822 2.022 2.604 0.615 1.433];
    [4.343 4.069 4.819 1.827 2.827 1.919 4.583 2.432 3.189 0.805 1.962];
    [0.189 0.830 1.131 0.256 1.446 0.225 0.379 1.167 0.766 0.799 0.273];
    [0.467 0.719 0.726 0.014 0.311 0.105 0.654 0.409 2.214 0.0002 2.648];
    [0.117 0.125 0.172 2.145 3.904 6.983 6.637 0.119 2.891 4.109 4.555];
    [9.809 6.137 8.136 4.495 4.318 6.997 8.825 5.987 8.236 19.86 5.442];
    [1.170 1.246 1.719 1.735 0.893 4.903 0.510 1.331 1.146 2.854 1.912];
    [4.000 1.648 0.0 0 0 0 0 0 0 0 0];
    [0.052 0.179 0.0 0 0 0 0 0 0 0 0];
    [18.56 27.97 0.0 0 0 0 0 0 0 0 0];
    [0.516 1.795 0.0 0 0 0 0 0 0 0 0]
    ]

    # select the right material, make sure it exists
  idx = findfirst(x-> x .== metal, metals)
  omega = 2*pi*c_const/lambda0
  # chose appropriate model
  if model == "LD" || model == "D"
      materialdata = LDdata[:,idx]
      wp = materialdata[1]
      f0 = materialdata[2]
      gamma0 = materialdata[3]

      fj = collect(materialdata[4:3:end])
      gammaj = collect(materialdata[5:3:end])
      wj = collect(materialdata[6:3:end])

      wp = (e_const/ħ_const).*wp
      wj = (e_const/ħ_const).*wj
      gamma0 = (e_const/ħ_const).*gamma0
      gammaj = (e_const/ħ_const).*gammaj

      # drude model contributions
      epsilon_D = 1 - (f0*wp^2)/(omega*(omega + 1.0im*gamma0))

      if model == "D"
          epsilon = epsilon_D
      else
          epsilon = epsilon_D .+ sum((fj.*wp^2)./((wj.^2 .- omega^2) .- 1.0im*omega.*gammaj))
      end
      return epsilon
    # BB model
  elseif model == "BB"
      materialdata = BBdata[:,idx]
      # numerical constants from article
      wp = materialdata[1]
      f0 = materialdata[2]
      gamma0 = materialdata[3]
      fj = materialdata[4:4:end]
      gammaj = materialdata[5:4:end]
      wj = materialdata[6:4:end]
      sigmaj = materialdata[7:4:end]
      wp = (e_const/ħ_const).*wp
      gamma0 = (e_const/ħ_const).*gamma0
      gammaj .= (e_const/ħ_const).*gammaj
      wj = (e_const/ħ_const).*wj
      sigmaj = (e_const/ħ_const).*sigmaj
      omegap = sqrt(f0).*wp

      epsf = 1 - omegap^2.0/(omega*(omega+gamma0*1.0im))
      aj = sqrt.(omega^2 .+ 1.0im*omega.*gammaj)
      zplus = (aj+wj)./(sqrt(2).*sigmaj)
      zminus = (aj-wj)./(sqrt(2).*sigmaj)

      epsb = (1.0im*sqrt(pi).*fj.*wp^2/(2*sqrt(2).*aj.*sigmaj)).*(erfcx.(-1.0im*zplus)+erfcx.(-1.0im.*zminus))
      eps = epsf + sum(epsb)
      return eps
    end
end


function eps_AH(w, material)
    AH_name = ["LiF","LiCl","LiBr","LiI","NaF","NaCl","NaBr","NaI","KF","KCl","KBr","KI","RbF","RbCl","RbBr", "RbI"]
    eps_0   = [9.03,11.89, 13.28, 11.03, 5.10, 5.92, 6.29, 7.31, 5.51, 4.86, 4.92, 6.11, 6.50, 4.90, 4.87, 4.93]
    eps_inf = [1.93, 2.75, 3.16, 3.80, 1.74, 2.33, 2.60, 3.01, 1.85, 2.17, 2.36, 2.65, 1.93, 2.18, 2.34, 2.58]
    nu_TO   = [305, 207, 173, 142, 246.5, 164, 134, 116, 194, 142, 114, 102, 158, 116.5, 87.5, 75.5]
    gamma0   = [0.053, 0.053, 0.053, 0.053, 0.065, 0.055, 0.067, 0.180, 0.060, 0.064, 0.034, 0.048, 0.062, 0.076, 0.056, 0.058]
    cft = 0.12398*1e-3 # this converts cm^-1 to eV so that w must be given in eV
    ii = findfirst(x-> x == material, AH_name)
    eps0_ii = eps_0[ii]
    epsi_ii = eps_inf[ii]
    gam = gamma0[ii]
    wT = cft*nu_TO[ii]
    eps_w = epsi_ii + (eps0_ii-epsi_ii)/(1 - w^2/wT^2 - im*gam*(w/wT))
    return eps_w
end

function eps_phonon_test(w, material)
    # these are taken from  Kitel, chapter 14: Plasmons, Polaritons, Polarons (not sure if the values are correct, also the damping parameters not provided, I added them artificially)
    names = ["LiH","LiF","LiCl","LiBr","NaF","NaCl","NaBr","KF","KCl","KI","RbF","RbI","CsCl","CsI","TlCl","TlBr",
             "AgCl","AgBr","MgO","GaP","GaAs","GaSb","InP","InAs","InSb","SiC","C","Si","Ge"]
    eps_0   = [12.9, 8.9, 12, 13.2, 5.1, 5.9, 6.4, 5.5, 4.85, 5.1, 6.5, 5.5, 7.2, 5.65, 31.9, 29.8, 12.3, 13.1,
              9.8, 10.7, 13.13, 15.69, 12.37, 14.55, 17.88, 9.6, 5.5, 11.7, 15.8]
    eps_inf = [3.6, 1.9, 2.7, 3.2, 1.7, 2.25, 2.6, 1.5, 2.1, 2.7, 1.9, 2.6, 2.6, 3, 5.1, 5.4, 4, 4.6, 2.95,
               8.5, 10.9, 14.4, 9.6, 12.3, 15.6, 6.7, 5.5, 11.7, 15.8]
    w_TO   = [0.0724091329, 0.038179361, 0.0236975344, 0.0197479453, 0.029621918, 0.0204062102, 0.0164566211,
              0.0236975344, 0.0177731508, 0.012507032, 0.0190896805, 0.0092157078, 0.012507032, 0.0078991781,
              0.0078991781, 0.0053319452, 0.012507032, 0.0098739727, 0.0493698633, 0.0454202743, 0.0335715071,
              0.0283053883, 0.0375210961, 0.0269888586, 0.0230392696, 0.0980814618, 0.165224476, 0.0651682196, 0.0375210961]

    gamma0   = [0.0586, 0.0669, 0.0588, 0.0622, 0.0533, 0.0557, 0.0506, 0.0531, 0.0697, 0.0595, 0.0586, 0.0669, 0.0588, 0.0622, 0.0533,
                0.0557, 0.0506, 0.0531, 0.0697, 0.0595, 0.0586, 0.0669, 0.0588, 0.0622, 0.0533, 0.0557, 0.0506, 0.0531, 0.0697]

    ii = findfirst(x -> x == material, names)
    eps0_ii = eps_0[ii]
    epsi_ii = eps_inf[ii]
    γ = gamma0[ii]
    wT = w_TO[ii]
    eps_w = epsi_ii + (eps0_ii-epsi_ii)/(1 - w^2/wT^2 - im*γ*(w/wT))
    return eps_w
end

function eps_VW(w, name)
    #from Martin Weismann1 and Nicolae C. Panoiu
    name_list = ["WS2", "WSe2", "MoS2", "MoSe2"]

    E = [2.009 2.204 2.198 2.407 2.4 2.595 2.644 2.831 3.056 3.577 5.078 5.594 ;
          1.654 2.426 2.062 2.887 2.2 2.6 3.8 5 NaN NaN NaN NaN;
          1.866 2.005 2.862 2.275 3.745 NaN NaN NaN NaN NaN NaN NaN;
          1.548 1.751 2.151 2.609 3.959 NaN NaN NaN NaN NaN NaN NaN ]

    f = [1.928 0.197 0.176 0.142 2.98 0.54 0.05 12.6 8.765 29.99 49.99 79.99;
          0.557 5.683 1.036 16.11 1.5 1.5 70 80 NaN NaN NaN NaN;
          0.752 1.883 36.89 10 100 NaN NaN NaN NaN NaN NaN NaN;
          0.648 1.302 4.621 37.4 121.4 NaN NaN NaN NaN NaN NaN NaN ]

    γ = [0.032 0.25 0.161 0.112 0.167 0.213 0.171 0.266 0.24 1.196 1.9 2.51;
          0.036 0.243 0.115 0.344 0.3 0.3 0.7 0.7 NaN NaN NaN NaN;
          0.045 0.097 0.383 1.0 0.533 NaN NaN NaN NaN NaN NaN NaN;
          0.043 0.097 0.537 0.582 0.896 NaN NaN NaN NaN NaN NaN NaN; ]
    ii = findfirst(name_list, name)
    eps = 1.0
    for j = 1:12
        if !isnan(E[ii, j])
            eps += f[ii, j]/(E[ii, j]^2 - w^2- 1.0im*w*γ[ii, j])
        end
    end
    return eps
end


function eps_pJC(w, name)
    # reference :
    # "Optimizing the Drude-Lorentz model for material permittivity:
    #    Method, program, and examples for gold, silver, and copper"
    #    H. S. Sehmi, W. Langbein, and E. A. Muljarov
    names = ["Ag", "Cu", "Au"]
    einf = [0.77259, 12.294, 0.83409]
    gamD = [0.02228, 0.07044, 0.02334]
    sigD = [3751.4, 1137.9, 3134.5]

    Omega_p = [3.9173 2.1508 2.6905; 3.988 4.6366 2.8772; 4.0746 4.9297 3.7911; 4.6198 8.8317 4.8532]
    Omega_pp = [-0.06084 -0.23449 -0.16645; -0.04605 -0.68811 -0.44473; -0.63141 -4.6932 -0.81981; -2.8279 -0.2679 -13.891]
    sig_p = [0.09267 0.95283 -0.01743; -0.0015342 0.97953 1.0349; 1.4911 -61.583 1.2274; 4.2843 -12.186 9.85]
    sig_pp = [0.01042 -0.12983 0.3059; -0.062233 0.48395 1.2919; 0.40655 35.021 2.5605; 4.2181 5.1474 -37.614]

    ii = findfirst(x-> x .== name, names)

    epsc = einf[ii] .- gamD[ii]*sigD[ii]./(w.*(w .+ 1.0im*gamD[ii]))
    for k = 1:4
        epsc = epsc .+ 1.0im*((sig_p[k, ii] + 1.0im*sig_pp[k, ii])./(w .- (Omega_p[k, ii] .+ 1.0im*Omega_pp[k, ii]))
                      +(sig_p[k, ii] - 1.0im*sig_pp[k, ii])./(w .+ (Omega_p[k, ii] .- 1.0im*Omega_pp[k, ii])))
    end
    return epsc
end

function eps_JC_vecfit(w, metal)
    metals = ["Au", "Ag1", "Ag2", "Cu", "Cr", "Al"]

    e = [1.0, 1.0, 1.0, 1.7665, 2.1508, 1.0]
    d = [9.6320E+2,  4.6904E+4, 2.5926E+3, 7.4394E+2, 1.7292E+2,  1.3186E+3]
    c = [ -9.6287e2  (7.2435e-1-1.0im*1.4189)  (9.7416+1.0im*3.4040e-2) (-7.7534e-1- 1.0im*4.2546e-1);
         -4.6896E+4 (-1.0859 + 1.0im*9.0177) (2.8247E-1 + 1.0im*3.3354E-1) (3.6945 - 1.0im*7.1655E-1);
         -2.5956E+3 1.1961E+1 ( 1.0284E-1 + 1.0im*3.999E-1) (3.1782 - 1.0im*5.5464E-1);
         -7.4532E+2 (1.0694 + 1.0im*1.5767E-1) (2.6922 - 1.0im*9.6133) (6.6689E-1 - 1.0im*1.4793);
         -1.6340E+2 (-7.1487 - 1.0im*8.8557E+1) (-1.8454 + 1.0im*1.1286E-1) (2.2579 - 1.0im*1.0288)
         -1.3541E+3 (1.3445E+1 - 1.0im*3.7366E+1) (4.1114 - 1.0im*5.8421) 0.0]


    p = [-8.0129E-2  (-4.7259E-1 + 1.0im*2.7018) (-1.8003 + 1.0im*2.7991) ( -7.4424E-1 + 1.0im*4.4474);
         -1.7367E-3 (-9.6710E-1 + 1.0im*6.2681E-1) (-2.5138E-1 + 1.0im*3.9208) (-1.0766 + 1.0im*3.9631);
         -3.4089E-2 -2.4860 (-2.5434E-1 + 1.0im*3.8737) (-8.9100E-1 + 1.0im*3.9425);
         -1.0616E-1 (-2.5740E-1 + 1.0im*2.1233) (-2.1296 + 1.0im*2.6916) (-8.3075E-1 + 1.0im*4.8844);
         -1.4642E-1 (-1.8013 + 1.0im*1.1172) (-3.3664E-1 + 1.0im*2.6614) (-1.1450 + 1.0im*5.4642);
         -1.4661E-1 (-9.2093E-1 + 1.0im*7.4420E-1) (-2.0507E-1 + 1.0im*1.4871) 0.0]

    ii = findfirst(x-> x .== metal, metals)
    s = 1.0im*w
    epsc = e[ii] .+ d[ii]./s
    for j = 1:size(p, 2)
        epsc = epsc .+ c[ii, j]./(s .- p[ii, j])
    end
    return conj(epsc)
end


function CPDOS(w, eps_inf, wp, wg, A1, A2, w01, w02, gamp, gam1, gam2)
    # Reference : A n Analytic Model for the Dielectric Function of Au, Ag, and their Alloys
    # by D. Rioux, S. Vallières et. al.
    eps_drude = eps_inf - wp^2/(w*(w + 1.0im*gamp))
    im_eps = -0.5*A2*log(1 - ((w + 1.0im*gam2)/w02)^2)/(w + 1.0im*gam2)^2
    re_eps = -0.5*sqrt(wg - w01)*log(1 - ((w + 1.0im*gam1)/w01)^2)/(w + 1.0im*gam1)^2
    re_eps = re_eps + 2*sqrt(wg)*atanh(sqrt((wg - w01)/wg))/(w + 1.0im*gam1)^2
    re_eps = re_eps - sqrt(w + 1.0im*gam1 - wg)*atan(sqrt((wg - w01)/(w + 1.0im*gam1 - wg)))/(w + 1.0im*gam2)^2
    re_Esp = re_eps - sqrt(w + 1.0im*gam1 + wg)*atanh(sqrt((wg - w01)/(w + 1.0im*gam1 + wg)))/(w + 1.0im*gam1)^2

    return eps_drude + re_eps + 1.0im*im_eps
end

function AgAu_eps_CPDOS(w, name)
    if name == "Au"
        return CPDOS(w, 2.2715, 8.9234, 2.6652, 73.251, 40.007, 2.3957, 3.5362, 0.042389, 0.1788, 0.35467)
    elseif name == "Ag"
        return CPDOS(w, 1.738, 8.5546, 4.0575, 51.217, 30.770, 3.9260, 4.1655, 0.022427, 0.017723, 0.18819)
    end
end



function nZEP(lam)
    x = 10*lam
    n0 = 1.541093
    n1 = 4.113002e5
    n2 = 4.070357e12
    n = n0 .+ n1./x.^2 + n2./x.^4
    return n
end


function epshBN_JGA(w)
  # BNx  parallel to film
  wl = 1.23984193./w
  ww = 10000 ./wl               # ww in cm-1

  exx = 4.95 .+ 1.23e5./(767^2 .- ww.*(ww .+ 1.0im*35)) .+ 3.49e6./(1367^2 .- ww.*(ww .+ 1.0im*29))
  #BNz (=BN)     perp to film
  ezz =  4.1 .+ 3.26e5./(783^2 .- ww.*(ww .+ 1.0im*8)) .+ 1.04e6./(1510^2 .- ww.*(ww .+ 1.0im*80))
  return [exx exx ezz]
end

function epshBN_Cai(w; id="clean")
    @assert id == "clean" || id == "damaged"
    # w is in eV
    ww = w*1e3 # ww is in meV to match the units of parameters
    if id == "clean"
        epsxx = 4.87 .+ 1.83*170.1^2 ./(170.1^2 .- ww.*(ww .+ 1.0im*0.87))
        epszz = 2.95 .+ 0.61*92.5^2 ./(92.5^2 .- ww.*(ww .+ 1.0im*0.25))
    else
        epsxx = 4.87 .+ 1.83*170.1^2 ./(170.1^2 .- ww.*(ww .+ 1.0im*6.5))
        epszz = 2.95 .+ 0.61*92.5^2 ./(92.5^2 .- ww.*(ww .+ 1.0im*1.9))
    end
    return [epsxx epsxx epszz]
end


function eps_hBN_Sandra(w)
    w = w*1000 # convert ev to meV
    epszz = 4.1 .+ 70.8.^2 ./(97.1^2 .- w.*(w .+ 0.99im)) .+ 126^2 ./(187^2 .- w.*(w .+ 9.92im))
    epsxx = 4.95 .+ 232^2 ./(170^2 .- w.*(w .+ 3.6im)) .+ 43.5^2 ./(95.1^2 .- w.*(w .+ 3.4im))
    return [epsxx epsxx epszz]
end


function eps_hBN_Andrea(w)
    wl = 1.23984193./w
    ww = 10000 ./wl               # ww in cm-1
    exx = 4.9.*(1.0 .+ (1614^2 - 1360^2)./(1360^2 .- ww.^2 .- 7.0im.*ww))
    #BNz (=BN)     perp to film
    ezz = 2.95.*(1.0 .+ (825^2 - 760^2)./(760^2 .- ww.^2 .- 2.0im.*ww))
    return [exx exx ezz]
end



function epsSi3N4(ω)
    # reference : Infrared dielectric properties of low-stress silicon nitride
    #             Giuseppe Cataldo et. al.
    # ω in eV
    ω = ω/0.00413566553853599  # makes the ω[THz]
    p = [7.582 0.0000 13.913 5.810 0.0001;
         6.754 0.3759 15.053 6.436 0.3427;
         6.601 0.0041 24.521 2.751 0.0006;
         5.430 0.1179 26.440 3.482 0.0002;
         4.601 0.2073 31.724 5.948 0.0080;
         4.562 0.0124 0.0000 0.000 0.0000]
    epsc = complex(p[6, 1], p[6, 2])     # ϵinf = eps[M+1], M = 5
    for j = 1:5
        Γ = p[j, 4]*exp(-p[j, 5]*((p[j, 3]^2 - ω^2)/(ω*p[j, 4]))^2)
        epsc += complex(p[j, 1] - p[j+1, 1], p[j, 2] - p[j+1, 2])*p[j, 3]^2/(p[j, 3]^2 - ω^2 - 1.0im*ω*Γ)
    end
    return epsc
end


function PCB_Molecule(ω)
    s = [158, 164, 49]
    Ω = [1450, 1505, 1478]
    Γ = [8.3, 13.4, 6.0]
    ϵ = 2.8 .+ 0.0im
    for i = 1:3
        ϵ += s[i]^2 ./(Ω[i]^2 .- ω.^2 .- 1.0im*Γ[i].*ω)
    end
    return ϵ
end

function PEI_eps(ω, d)
    @assert 10 <= d <= 75
    Ω1 = [1300 0.012 0.05; 1404 0.005 0.04; 1474 0.013 0.04;   # 10 nm
          1560 0.019 0.04; 1700 0.005 0.06]
    ε1 = 2.56 + 0.0im
    Ω2 = [1348.7 0.0157 0.0555; 1465 0.0126 0.0493; 1576 0.0136 0.0524]
    ε2 = 2.0717 + 0.0im
    ε = (ε2-ε1)/65*d + (15*ε1 - 2*ε2)/13

    for i = 1:3
        S =  (Ω2[i, 2]-Ω1[i, 2])/65*d + (15*Ω1[i, 2] - 2*Ω2[i, 2])/13
        Ω =  (Ω2[i, 1]-Ω1[i, 1])/65*d + (15*Ω1[i, 1] - 2*Ω2[i, 1])/13
        γ =  (Ω2[i, 3]-Ω1[i, 3])/65*d + (15*Ω1[i, 3] - 2*Ω2[i, 3])/13
        ε += S./(1 .- (ω./Ω).^2 .- 1.0im*γ.*(ω./Ω))
    end
    return ε
end


function interp_PEI_params(d)
    @assert 10 <= d <= 75
    Ω1 = [1300 0.012 0.05; 1404 0.005 0.04; 1474 0.013 0.04;   # 10 nm
          1560 0.019 0.04; 1700 0.005 0.06]
    ε1 = 2.56 + 0.0im
    Ω2 = [1348.7 0.0157 0.0555; 1465 0.0126 0.0493; 1576 0.0136 0.0524]
    ε2 = 2.0717 + 0.0im
    ε = (ε2-ε1)/65*d + (15*ε1 - 2*ε2)/13
    S = zeros(3)
    Ω = zeros(3)
    γ = zeros(3)
    for i = 1:3
        S[i] =  (Ω2[i, 2]-Ω1[i, 2])/65*d + (15*Ω1[i, 2] - 2*Ω2[i, 2])/13
        Ω[i] =  (Ω2[i, 1]-Ω1[i, 1])/65*d + (15*Ω1[i, 1] - 2*Ω2[i, 1])/13
        γ[i] =  (Ω2[i, 3]-Ω1[i, 3])/65*d + (15*Ω1[i, 3] - 2*Ω2[i, 3])/13
    end
    X = [Ω[1]; S[1]; γ[1]; Ω[2]; S[2]; γ[2]; Ω[3]; S[3]; γ[3]; real(ε)]
    return X
end





#=
ω = range(1000, 3000, length =500
d = range(10, 75, step = 0.1)
epsc = zeros(ComplexF64, 5000, length(d))
for (j, dj) = enumerate(d)
    for (i, ωi) = enumerate(ω)
        #epsc[i] = PCB_Molecule(ωi)
        epsc[i, j] = PEI_eps(ωi, dj)
    end
end

using PyPlot
pygui(true)

pcolormesh(ω, d, imag.(epsc)', cmap = "hot")

begin
    for j = 1:100:length(d)
        #plot(ω, real.(epsc[:, j]))
        plot(ω, imag.(epsc[:, j]))
    end
end

using DelimitedFiles
writedlm("C:/Users/vmkhitaryan/Documents/Github/jcalc/eps_data/Si3N4(nk-eV-Analytic).txt", [ω real.(nc) imag.(nc)])

=#
