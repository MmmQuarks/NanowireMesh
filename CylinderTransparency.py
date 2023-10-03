import numpy as np
from scipy.special import jv, jvp, hankel1, h1vp

def _series_approximation(func, nStart, tolerance):
	assert 0 < tolerance < 1, 'Tolerance must be between 0 and 1'

	# initializing the list of terms with the first term
	n = nStart
	terms = [func(n)]
	relativeChange = 1
	while relativeChange > tolerance or n < 50:
		n += 1
		newTerm = func(n)
		relativeChange = ( np.abs(newTerm) - np.abs(sum(terms)) ) / np.abs(sum(terms))
		terms.append(newTerm)
	print('iterations:', n)
	return sum(terms)

# a little syntactic sugar
def Jn(n, z): # Bessel function of first kind
	return jv(n, z)

def Jnp(n, z): # First deriv of Bessel function of first kind
	return jvp(n, z)

def Hn(n, z): # Hankel function of first kind
	return hankel1(n, z)

def Hnp(n, z): # First deriv of Hankel function of first kind
	return h1vp(n, z)

def extinction_coefficient(wavelength, length, mr, xr):
	# a terms are for C_s
	# b terms are for C_p
	real = np.real

	def _a_n(n, mr = mr, xr = xr): # this yields Real{a_n + b_n}
		aNumerator =  Jn(n, xr) * Jnp(n, mr * xr) - mr * Jnp(n, xr) * Jn(n, mr * xr) 
		aDenominator = Hn(n, xr) * Jnp(n, mr * xr) - mr * Hnp(n, xr) * Jn(n, mr * xr)
		a_n = aNumerator / aDenominator
		return real(a_n)

	def _b_n(n, mr = mr, xr = xr): # this yields Real{b_n}
		bNumerator =  mr * Jn(n, xr) * Jnp(n, mr * xr) - Jnp(n, xr) * Jn(n, mr * xr) 
		bDenominator = mr * Hn(n, xr) * Jnp(n, mr * xr) -  Hnp(n, xr) * Jn(n, mr * xr)
		b_n = bNumerator / bDenominator
		return real(b_n)

	def _nthTerm(n, mr = mr, xr = xr): # this yields Real{a_n + b_n}
#		aNumerator =  Jn(n, xr) * Jnp(n, mr * xr) - mr * Jnp(n, xr) * Jn(n, mr * xr) 
#		aDenominator = Hn(n, xr) * Jnp(n, mr * xr) - mr * Hnp(n, xr) * Jn(n, mr * xr)
#		a_n = aNumerator / aDenominator
#
#		bNumerator =  mr * Jn(n, xr) * Jnp(n, mr * xr) - Jnp(n, xr) * Jn(n, mr * xr) 
#		bDenominator = mr * Hn(n, xr) * Jnp(n, mr * xr) -  Hnp(n, xr) * Jn(n, mr * xr)
#		b_n = bNumerator / bDenominator

		return _a_n(n, mr = mr, xr = xr) + _b_n(n, mr = mr, xr = xr)


	# the sum of Real{a_0 + b_0}
	a0b0 = _nthTerm(n = 0)
	a0 = _a_n(n = 0)
	b0 = _b_n(n = 0)

	# now calculating sum for Real(a_n + b_n) for 1 < n < infty
	anbn = _series_approximation(func = _nthTerm, nStart = 1, tolerance = 0.000001)
	an = _series_approximation(func = _a_n, nStart = 1, tolerance = 0.000001)
	bn = _series_approximation(func = _b_n, nStart = 1, tolerance = 0.000001)

	C_s = 2 * wavelength * length / np.pi * (a0 + 2 * an)
	C_p = 2 *  wavelength * length / np.pi * (b0 + 2 * bn)
	C_ext = wavelength * length /np.pi * (a0b0 + 2 * anbn)

	C_ext = 1/2 * (C_s + C_p)
	return dict(C_s = C_s, C_p = C_p, C_ext = C_ext) 

if __name__ == "__main__":
	mr = 0.055 + 3.32j # refract index of silver over vacuum
	wavelength = 550E-9 # wavelength in vacuum
	k = 2 * np.pi / wavelength # light wavenumber in vacuum
	R = (.15E-6) / 2 # cylinder radius
	xr = k * R # size parameter
	L = 10E-6
	ans = extinction_coefficient(wavelength = wavelength, length = L, mr = mr, xr = xr) 
	print(ans)

	print('calculating transparency of network at 50 amd')
	amd = 50 # mg / m^2
	silverDensity = 1.049E10 #(mg/m^3)
	wireMass = np.pi * R**2 * L * silverDensity # mass in mg
	print('wire mass =', wireMass)
	n_s = amd / wireMass
	print('n_s', n_s)
	T = np.exp(- n_s * ans['C_ext'])
	print(T)
