using SpecialFunctions


export besseljp, besselyp, hankelh1p, hankelh2p, besselkp, sp_besselj, sp_besseljp, sp_bessely, sp_besselyp
export sp_hankelh1, sp_hankelh1p, sp_hankelh2, sp_hankelh2p, RiccatiPsi, RiccatiPsip, RiccatiXi, RiccatiXip
export RiccatiChi, RiccatiChip

###############    Bessel function Derivatives (Probably can be done better with macros) ############################
function besseljp(n::Real, x::Number)
  if n == 0
    return -besselj(1,x)
  else
    besselj(n-1, x) - (convert(typeof(x), n)/x)*besselj(n, x)
  end
end

function besselyp(n::Real, x::Number)
  if n == 0
    return -bessley(1,x)
  else
    bessely(n-1, x) - (convert(typeof(x), n)/x)*bessely(n, x)
  end
end

function hankelh1p(n::Real, x::Number)
  if n== 0
    return -hankelh1(1, x)
  else
    hankelh1(n-1, x) - (convert(typeof(x), n)/x)*hankelh1(n, x)
  end
end

function hankelh2p(n::Real, x::Number)
  if n == 0
    return -hankelh2(1, x)
  else
    hankelh2(n-1, x) - (convert(typeof(x), n)/x)*hankelh2(n, x)
  end
end

function besselkp(n::Real, x::Number)
	return -(besselk(n-1, x) + (n/x)*besselk(n, x))
end


# Define the modified Bessel Function of the second kind with complex order & argument
function CBesselK(n::Number, z::Number)

	Knu = exp(lgamma(nu))  * (z/2)^(-nu)
	Kmu = exp(lgamma(-nu)) * (z/2)^nu

	bound = 20

	precision = 1e-10

	s = 0

	if (abs(z) < 10)
		for i = 0:bound
			var = i
			facto = factorial(i)

			r = 0.5 * Knu * ((z/2)^(2*var) * exp(-lgamma(1+var-nu)+lgamma(1-nu))/facto) + 0.5 * Kmu * ((z/2)^(2*var) * exp(-lgamma(1+var+nu)+lgamma(1+nu))/facto)
			s+=r

			if ((abs(real(r/s))<precision)&&(abs(imag(r/s))<precision))
				break
			end
		end
	else
		for i = 0:bound
			var = i
			facto = factorial(i)

			r = sqrt(pi/(2*z)) * exp(-z) * exp(lgamma(0.5+var+nu)+lgamma(0.5+var-nu)-lgamma(0.5+nu)-lgamma(0.5-nu)/facto) * (-1/2*z)^var
			s+=r

			if ((abs(real(r/s))<precision)&&(abs(imag(r/s))<precision))
				break
			end
		end
	end

	return s

end

################### Spherical bessel functions and their derivatives ################

function sp_besselj(n::Real, x::Number)
    return sqrt(0.5*pi/x)*besselj(n+0.5, x)
end

function sp_besseljp(n::Real, x::Number)
  return -sp_besselj(n+1,x) + (convert(typeof(x), n)/x)*sp_besselj(n,x)
end
################
function sp_bessely(n::Real, x::Number)
  return sqrt(0.5*pi/x)*bessely(n+0.5, x)
end
function sp_besselyp(n::Real, x::Number)
  return -sp_bessely(n+1,x) + (convert(typeof(x), n)/x)*sp_bessely(n,x)
end
#################

function sp_hankelh1(n::Real, x::Number)
    return sqrt(0.5*pi/x)*hankelh1(n+0.5,x)
end
function sp_hankelh1p(n::Real, x::Number)
  return -sp_hankelh1(n+1,x) + (convert(typeof(x), n)/x)*sp_hankelh1(n,x)
end
#################
function sp_hankelh2(n::Real, x::Number)
  return sqrt(0.5*pi/x)*hankelh2(n+0.5,x)
end
function sp_hankelh2p(n::Real, x::Number)
  return -sp_hankelh2(n+1,x) + (convert(typeof(x), n)/x)*sp_hankelh2(n,x)
end

###################### Riccati-Bessel ##################

function RiccatiPsi(n::Real, z::Number)  # naive way of calculating
  return sqrt(0.5*pi*z)*besselj(n+0.5, z)
end
function RiccatiPsip(n::Real, z::Number)
  return ((n+1)/z)*RiccatiPsi(n, z)-RiccatiPsi(n+1, z)
end
###############
function RiccatiXi(n::Real, z::Number)
  return sqrt(0.5*pi*z)*hankelh1(n+0.5, z)
end
function RiccatiXip(n::Real, z::Number)
  return ((n+1)/z)*RiccatiXi(n, z)- RiccatiXi(n+1, z)
end
########################
function RiccatiChi(n::Real, z::Number)
  return -sqrt(0.5*pi*z)*bessely(n+0.5, z)
end

function RiccatiChip(n::Real, z::Number)
  return -((n+1)/z)*RiccatiChi(n,z) + RiccatiChi(n+1,z)
end
