from main import *

n = 1
rhob = n*1.22*PROTON
rhodm = rhob * (0.315-0.048)/0.048
Tb = 2e2
Tdm = 10
uth = uthf(PROTON, 0.3*GeV_to_mass, Tb, Tdm)/1e5
lv = np.logspace(np.log10(uth)-1, np.log10(uth)+2, 100)
lQ = [Q_IDMB(rhodm, v*1e5, Tb, Tdm)*BOL for v in lv]

x1, x2 = uth/10, uth*100
y1, y2 = np.min(lQ)*1.05, np.max(lQ)*1.5
plt.figure()
plt.plot(lv, lQ)
plt.plot([x1, x2], [0, 0], 'k', lw = 0.5)
plt.plot([uth, uth], [y1, y2], 'k--', label=r'$u_{th}$')
plt.legend()
plt.xscale('log')
plt.xlim(x1, x2)
plt.ylim(y1, y2)
plt.xlabel(r'$v_{b\chi}\ [\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$(\gamma-1)\dot{Q}_{b}\ [\mathrm{erg\ s^{-1}}]$')
plt.text(uth*1.1, np.min(lQ)*0.9, r'$T_{b}=200\ \mathrm{K}$, $T_{\chi}=10\ \mathrm{K}$, $n=1\ \mathrm{cm^{-3}}$')
plt.tight_layout()
plt.savefig('Qb_v.pdf')
#plt.show()

lv = [0, 15, 30, 60]
lls = ['-', '--', '-.', ':']
lla = [r'$v_{b\chi,0}=0$', r'$v_{b\chi,0}=0.5\sigma_{rms}$', r'$v_{b\chi,0}=1\sigma_{rms}$', r'$v_{b\chi,0}=2\sigma_{rms}$']

rep = '100-sigma_ns/'

yup = 2e8
plt.figure()
for v0, ls, lab in zip(lv, lls, lla):
	lm_, lz_, lxh2_, lxhd_, lxe_, lTb_, lvr_ = np.array(retxt(rep+'Mthz_BDMS_'+str(v0)+'.txt',7,0,0))
	plt.plot(lz_, lm_, label=lab+', BDMS', ls = ls, lw = 2, color='r')
for v0, ls, lab in zip(lv, lls, lla):
	lm, lz, lxh2, lxhd, lxe, lTb, lvr = np.array(retxt(rep+'Mthz_CDM_'+str(v0)+'.txt',7,0,0))	
	plt.plot(lz, lm, label=lab+', CDM', ls = ls, lw = 1, color='b')
plt.plot(lz, Mdown(lz), 'k-', label='Trenti & Stiavelli (2009)', lw = 0.5)
plt.fill_between([15,20],[1e4,1e4],[yup,yup], facecolor='gray', label='EDGES')
plt.legend()
plt.xlabel(r'$z_{vir}$')
plt.ylabel(r'$M_{\mathrm{th}}\ [M_{\odot}]$')
#plt.xscale('log')
plt.yscale('log')
plt.xlim(15, 100)
plt.ylim(1e4, yup)
plt.tight_layout()
plt.savefig('Mth_z.pdf')
plt.close()

m = 1e6
zvir = 20
#rep0 = 'example/'
for v0 in lv:
	d0 = readd(m, zvir, v0, mode = 0)
	d1 = readd(m, zvir, v0, mode = 1)
	plt.figure()
	plt.plot(d0['t'], d0['Tb'], label=r'$T_{b}$, CDM', color='b')
	plt.plot(d1['t'], d1['Tb'], '--', label=r'$T_{b}$, BDMS', color='r')
	#plt.plot(d0['t'], d0['Tdm'], '-.', label=r'$T_{\chi}$, CDM')
	plt.plot(d1['t'], d1['Tdm'], '-.', label=r'$T_{\chi}$, BDMS', color='orange')
	plt.plot(d1['t'], (d1['z']+1)*2.726, ':', label=r'$T_{cmb}$', color='g')
	plt.xlabel(r'$t\ [\mathrm{Myr}]$')
	plt.ylabel(r'$T\ [\mathrm{K}]$')
	plt.legend()
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(np.min(d0['t']), np.max(d0['t']))
	plt.ylim(1, 3e3)
	plt.tight_layout()
	plt.savefig('Example_T_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')
	plt.close()

#lv = [0 15, 30, 60]
#for v0 in lv:
#	d0 = readd(m, zvir, v0, mode = 0)
#	d1 = readd(m, zvir, v0, mode = 1)
	if v0==0:
		continue
	vIGM = [vbdm_z(z, v0)/1e5 for z in d0['z']]
	luth = [uthf(PROTON, 0.3*GeV_to_mass, T, Td)/1e5 for T, Td in zip(d1['Tb'], d1['Tdm'])]
	plt.figure()
	plt.plot(d0['t'], d0['v']/1e5, label='CDM', color='b')
	plt.plot(d1['t'], d1['v']/1e5, '--', label='BDMS', color='r')
	plt.plot(d1['t'], luth, '-.', label=r'$u_{th}$', color='orange')
	plt.plot(d0['t'], vIGM, ':', label='IGM, CDM', color='g')
	plt.xlabel(r'$t\ [\mathrm{Myr}]$')
	plt.ylabel(r'$v_{b\chi}\ [\mathrm{km\ s^{-1}}]$')
	plt.legend()
	plt.ylim(np.min(vIGM)/10, np.max(vIGM)*10)
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig('Example_vbDM_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')
	plt.close()

lv = np.logspace(0, 3, 31)
vd, vu = np.min(lv), np.max(lv)
lz = retxt(rep+'zbase.txt',1,0,0)[0]
nz = len(lz)

lm = np.array(retxt(rep+'Mth_v.txt',nz,0,0))
lxh2 = np.array(retxt(rep+'xh2_v.txt',nz,0,0))
lvr = np.array(retxt(rep+'vr_v.txt',nz,0,0))
lm_ = np.array(retxt(rep+'Mth_v0.txt',nz,0,0))
lxh2_ = np.array(retxt(rep+'xh2_v0.txt',nz,0,0))
lvr_ = np.array(retxt(rep+'vr_v0.txt',nz,0,0))

y1, y2 = 1e4, 3e9
plt.figure()
a = [plt.plot(lv, lm[i], label=r'$z_{vir}='+str(int(lz[i]*100)/100)+'$, BDMS', ls = lls[i], color='r', lw = 2) for i in range(nz)]
a = [plt.plot(lv, lm_[i], label=r'$z_{vir}='+str(int(lz[i]*100)/100)+'$, CDM',ls=lls[i], color='b', lw=1) for i in range(nz)]
plt.plot([30,30],[y1,y2], 'k', lw=0.5)
plt.plot([60,60],[y1,y2], 'k', lw=0.5)
plt.plot([90,90],[y1,y2], 'k', lw=0.5)
plt.legend(loc=2)
plt.xlabel(r'$v_{b\chi,0}\ [\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$M_{\mathrm{th}}\ [M_{\odot}]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(vd, vu)
plt.ylim(y1, y2)
plt.tight_layout()
plt.savefig('Mth_v.pdf')
plt.close()

y1, y2 = 4e-6, 7e-4
plt.figure()
a = [plt.plot(lv, lxh2[i], label=r'$z_{vir}='+str(int(lz[i]*100)/100)+'$, BDMS', ls = lls[i], color='r', lw = 2) for i in range(nz)]
a = [plt.plot(lv, lxh2_[i], label=r'$z_{vir}='+str(int(lz[i]*100)/100)+'$, CDM',color='b',ls=lls[i],lw=1) for i in range(nz)]
plt.plot([30,30],[y1,y2], 'k', lw=0.5)
plt.plot([60,60],[y1,y2], 'k', lw=0.5)
plt.plot([90,90],[y1,y2], 'k', lw=0.5)
plt.legend(loc=3)
plt.xlabel(r'$v_{b\chi,0}\ [\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$x_{\mathrm{H_{2}}}$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(vd, vu)
plt.ylim(y1, y2)
plt.tight_layout()
plt.savefig('XH2_v.pdf')
plt.close()

y1, y2 = 1e-1, 3e3
plt.figure()
a = [plt.plot(lv, lvr[i]/1e5, label=r'$z_{vir}='+str(int(lz[i]*100)/100)+'$, BDMS', ls = lls[i], color='r', lw=2) for i in range(nz)]
	#a = [plt.plot(lv, [vbdm_z(zr[i], v)/1e5  for v in lv], label=r'$z_{vir}='+str(int(lz[i]*100)/100)+'$, CDM',color='k',ls=lls[i],lw=0.5) for i in range(nz)]
a = [plt.plot(lv, lvr_[i]/1e5, label=r'$z_{vir}='+str(int(lz[i]*100)/100)+'$, CDM',color='b',ls=lls[i],lw=1) for i in range(nz)]
plt.plot([30,30],[y1,y2], 'k', lw=0.5)
plt.plot([60,60],[y1,y2], 'k', lw=0.5)
plt.plot([90,90],[y1,y2], 'k', lw=0.5)
plt.legend()
plt.xlabel(r'$v_{b\chi,0}\ [\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$v_{b\chi,V}\ [\mathrm{km\ s^{-1}}]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(vd, vu)
plt.ylim(y1, y2)
plt.tight_layout()
plt.savefig('Vr_v.pdf')
plt.close()
