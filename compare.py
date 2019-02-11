from main import *

rep1 = '100-sigma_test/'
rep2 = '100-sigma_test0/'

nbin = 32
v0 = 24

X = np.array(retxt(rep1+'X_'+str(v0)+'.txt',nbin,0,0))
Y = np.array(retxt(rep1+'Y_'+str(v0)+'.txt',nbin,0,0))
Mh1 = np.array(retxt(rep1+'Mh_'+str(v0)+'.txt',nbin,0,0))
Mh2 = np.array(retxt(rep2+'Mh_'+str(v0)+'.txt',nbin,0,0))
dis = Mh2/Mh1 #np.abs(Mh2-Mh1)/Mh2

plt.figure()
ctf = plt.contourf(X, Y, dis, np.linspace(0, 2, nbin*2))
for c in ctf.collections:
	c.set_edgecolor('face')
cb = plt.colorbar()
cb.set_label(r'$\delta M_{\mathrm{th}}/M_{\mathrm{th}}$',size=12)
plt.plot([0.3], [8e-20], '*', color='purple')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$m_{\chi}c^{2}\ [\mathrm{GeV}]$')
plt.ylabel(r'$\sigma_{1}\ [\mathrm{cm^{2}}]$')
plt.tight_layout()
plt.savefig('Mth_com_v0_'+str(v0)+'.pdf')
plt.close()
