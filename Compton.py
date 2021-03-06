from main import *
from scipy.optimize import curve_fit

d0 = evolve(1e6, 1, mode = 0, dmax = delta0, v0 = 0, fac = 0)
lz = d0['z']
lxe = d0['X'][5][lz>20]
lz = lz[lz>20]
lg = np.array([GammaC(z) for z in lz])

def f(x, a, b):
	return a*x + b

fit = curve_fit(f, 4*np.log10(1+lz), np.log10(lg/lxe))
print(fit)

refa = 5.41e-36/BOL

plt.plot(4*np.log10(1+lz), np.log10(lg/lxe), label='Reference')
plt.plot(4*np.log10(1+lz), 4*np.log10(1+lz)*fit[0][0] + fit[0][1], '--', label='Fitting')
plt.plot(4*np.log10(1+lz), np.log10(refa*(1+lz)**4), '-.', label='Analytical')
plt.ylabel(r'$\log(\Gamma_{C}/x_{\mathrm{e}}\ [\mathrm{s^{-1}}])$')
plt.xlabel(r'$4\log(1+z)$')
plt.tight_layout()
plt.savefig('comption.pdf')

totxt('compton_ref.txt',fit,0,0,0)
