from txt import *
from scipy.interpolate import interp1d
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def JLW(rho, z, eta = 1e4):
	return 2*(eta/1e4)*(rho*1e2)*((1+z)/(10))**3

def Mratio(J21):
	return (1+6.96*(4*np.pi*J21)**0.47)

repref = 'boyuan_sfrd/'
refIII = np.array(retxt(repref+'popIII_sfr_s2_skx.txt',4,0,0))
refII = np.array(retxt(repref+'popII_sfr_s2_skx.txt',4,0,0))

lz1 = 1/refIII[0]-1
lz2 = 1/refII[0] -1

lJ1 = JLW(refIII[3], lz1)
lJ2 = JLW(refII[3], lz2, 2e3)

J_z1 = interp1d(lz1, lJ1)
J_z2 = interp1d(lz2, lJ2)

lJ0 = J_z1(lz1) + J_z2(lz1)

x1, x2 = 7, 27
down1, up1 = 1e-3, 2.
plt.figure()
plt.plot(lz1[lJ1>0], lJ1[lJ1>0], label='Pop III')
plt.plot(lz2[lJ2>0], lJ2[lJ2>0], '--', label='Pop II')
plt.plot(lz1[lJ0>0], lJ0[lJ0>0], '-.', label='Total (Jaacks et al. 2018)')
plt.plot([x1, x2], [0.1, 0.1], 'k:', label='Minihalo suppression')
plt.fill_between([15, 20],[up1, up1],[down1, down1],label='EDGES',facecolor='gray', alpha=0.5)
plt.fill_between([16, 19],[up1, up1],[down1, down1],facecolor='gray', alpha=0.5)
plt.xlabel(r'$z$')
plt.ylabel(r'$J_{\mathrm{LW},21}\equiv J_{\mathrm{LW}}/(10^{-21}\mathrm{erg\ s^{-1}\ cm^{-2}\ Hz^{-1}\ sr^{-1}})$')
plt.yscale('log')
plt.xlim(x1, x2)
plt.ylim(down1, up1)
plt.legend()
plt.tight_layout()
plt.savefig('J21_z.pdf')
#plt.show()
plt.close()

down2, up2 = 0, 30
plt.figure()
plt.plot(lz1[lz1<30], Mratio(lJ0)[lz1<30], label=r'$M_{\mathrm{th}}=M_{\mathrm{th},0}\left[1+6.96(4\pi J_{\mathrm{LW},21})^{0.47}\right]$'+'\n(Fialkov 2014)')
plt.fill_between([15, 20],[up2, up2],[down2, down2],label='EDGES',facecolor='gray', alpha=0.5)
plt.fill_between([16, 19],[up2, up2],[down2, down2],facecolor='gray', alpha=0.5)
plt.xlabel(r'$z$')
plt.ylabel(r'$M_{\mathrm{th}}/M_{\mathrm{th},0}$')
plt.xlim(x1, x2)
plt.ylim(down2, up2)
#plt.yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('Mrat_z.pdf')
plt.close()
