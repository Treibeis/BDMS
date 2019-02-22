from main import *
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

llogm = np.array([-4, -3, -2, -1.5, -1, -0.75, -0.5, -0.25, 0, 0.1, 0.2, 0.3, 0.43, 0.5, 0.6])
llogs = np.array([-0.45, -0.445, -0.44, -0.41, -0.35, -0.25, -0.15, 0, 0.2, 0.3, 0.5, 0.7, 1., 1.25, 4.])-20.

n = 1
rhob = n*1.22*PROTON
rhodm = rhob * (0.315-0.048)/0.048
Tb = 2e2
Tdm = 10
uth = uthf(PROTON, 0.3*GeV_to_mass, Tb, Tdm)/1e5
lv = np.logspace(np.log10(uth)-1, np.log10(uth)+2, 100)
lQ = [Q_IDMB(rhodm, v*1e5, Tb, Tdm)*BOL/(5/3-1) for v in lv]

x1, x2 = uth/10, uth*100
y1, y2 = np.min(lQ)*1.05, np.max(lQ)*1.5
plt.figure()
plt.plot(lv, lQ, 'g-')
plt.plot([x1, x2], [0, 0], 'k--', lw = 0.5)
plt.plot([uth, uth], [y1, y2], 'k-.', label=r'$u_{\mathrm{th}}$')
plt.fill_between([x1, x2], [0, 0], [y1, y1], facecolor = 'b', alpha=0.3, label='Cooling')
plt.fill_between([x1, x2], [0, 0], [y2, y2], facecolor = 'r', alpha=0.3, label='Heating')
plt.legend(loc=2)
plt.xscale('log')
plt.xlim(x1, x2)
plt.ylim(y1, y2)
plt.xlabel(r'$v_{b\chi}\ [\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$\dot{Q}_{b}\ [\mathrm{erg\ s^{-1}}]$')
plt.text(uth*1.1, np.min(lQ)*0.9, r'$T_{b}=200\ \mathrm{K}$, $T_{\chi}=10\ \mathrm{K}$, $n=1\ \mathrm{cm^{-3}}$')
plt.tight_layout()
plt.savefig('Qb_v.pdf')
plt.close()

lv = [0, 24, 45, 60]
lls = ['-', '--', '-.', ':']
lla = [r'$v_{b\chi,0}=0$', r'$v_{b\chi,0}=0.8\sigma_{\mathrm{rms}}$', r'$v_{b\chi,0}=1.5\sigma_{\mathrm{rms}}$', r'$v_{b\chi,0}=2\sigma_{\mathrm{rms}}$']

#rep = '100-sigma_ns/'
rep = '100-sigma_test0/'

#y0, yup = 5e4, 2e8
y0, yup = 5e4, 1e8
x1, x2 = 15, 60
lzt = np.linspace(x1, x2, 10)
lt = [TZ(z)/(1e6*YR) for z in lzt]
#y1, y2, y3 = 1e6*(1.63+0.48)/2, 1e6*(1.21+3.93)/2, 1e6*(11.27+6.08)/2
y1, y2, y3 = (1.63*0.48)**0.5*1e6, 1e6*(1.21*3.93)**0.5, 1e6*(11.27*6.08)**0.5
fig = plt.figure(figsize=(12,5.5))
ax1 = fig.add_subplot(121)
for v0, ls, lab in zip(lv[:2], lls[:2], lla[:2]):
	lm_, lz_, lxh2_, lxhd_, lxe_, lTb_, lvr_ = np.array(retxt(rep+'Mthz_BDMS_'+str(v0)+'.txt',7,0,0))
	plt.plot(lz_, lm_, label=lab+', BDMS', ls = ls, lw = 2, color='r')
for v0, ls, lab in zip(lv[:2], lls[:2], lla[:2]):
	lm, lz, lxh2, lxhd, lxe, lTb, lvr = np.array(retxt(rep+'Mthz_CDM_'+str(v0)+'.txt',7,0,0))	
	plt.plot(lz, lm, label=lab+', CDM', ls = ls, lw = 1, color='b')
plt.plot([x1, 24],[y1, y1], 'g-', label=lla[0]+', S18', lw = 3.5, alpha=0.8)
plt.plot([x1, 20],[y2, y2], 'g--', label=r'$v_{b\chi,0}=1\sigma_{\mathrm{rms}}$'+', S18', lw = 3.5, alpha=0.8)
plt.plot(lz, M_vcir(lz, Vcool(lz, 0)), 'k-', label=lla[0]+', F12', lw = 0.5)
plt.plot(lz, M_vcir(lz, Vcool(lz, lv[1])), 'k--', label=lla[1]+', F12', lw = 0.5)
plt.fill_between([15,20],[1e4,1e4],[yup,yup], facecolor='gray', label='EDGES')
plt.legend()
plt.xlabel(r'$z_{\mathrm{vir}}$')
plt.ylabel(r'$\tilde{M}_{\mathrm{th}}\ [M_{\odot}]$')
#plt.xscale('log')
plt.yscale('log')
plt.xlim(x1, x2)
plt.ylim(y0, yup)
ax2 = ax1.twiny()
ax2.set_xticks(lzt)
ax2.set_xticklabels([str(int(t*10)/10) for t in lt],size=14)
ax2.set_xlabel(r'$t\ [\mathrm{Myr}]$')
plt.xlim(x1, x2)

y0, yup = 1e5, 3e8
ax1 = fig.add_subplot(122)
for v0, ls, lab in zip(lv[2:], lls[2:], lla[2:]):
	lm_, lz_, lxh2_, lxhd_, lxe_, lTb_, lvr_ = np.array(retxt(rep+'Mthz_BDMS_'+str(v0)+'.txt',7,0,0))
	plt.plot(lz_, lm_, label=lab+', BDMS', ls = ls, lw = 2, color='r')
for v0, ls, lab in zip(lv[2:], lls[2:], lla[2:]):
	lm, lz, lxh2, lxhd, lxe, lTb, lvr = np.array(retxt(rep+'Mthz_CDM_'+str(v0)+'.txt',7,0,0))	
	plt.plot(lz, lm, label=lab+', CDM', ls = ls, lw = 1, color='b')
plt.plot([x1, 18],[y3, y3], 'g:', label=lla[3]+', S18', lw = 3.5, alpha=0.8)
plt.plot(lz, M_vcir(lz, Vcool(lz, lv[2])), 'k-.', label=lla[2]+', F12', lw = 0.5)
plt.plot(lz, M_vcir(lz, Vcool(lz, lv[3])), 'k:', label=lla[3]+', F12', lw = 0.5)
plt.fill_between([15,20],[1e4,1e4],[yup,yup], facecolor='gray')#, label='EDGES')
plt.legend()
plt.xlabel(r'$z_{\mathrm{vir}}$')
plt.ylabel(r'$\tilde{M}_{\mathrm{th}}\ [M_{\odot}]$')
plt.yscale('log')
plt.xlim(x1, x2)
plt.ylim(y0, yup)
ax2 = ax1.twiny()
ax2.set_xticks(lzt)
ax2.set_xticklabels([str(int(t*10)/10) for t in lt],size=14)
ax2.set_xlabel(r'$t\ [\mathrm{Myr}]$')
plt.xlim(x1, x2)
plt.tight_layout()
plt.savefig('Mth_z.pdf')
plt.close()

cons = 1.7e-17*GeV_to_mass * (SPEEDOFLIGHT/1e5)**4

m = 1e6
zvir = 20
nbin = 64
drep = 'data0/'
#rep0 = 'example/'
for v0 in lv:
	X = np.array(retxt(rep+'X_'+str(v0)+'.txt',nbin,0,0))
	Y = np.array(retxt(rep+'Y_'+str(v0)+'.txt',nbin,0,0))
	Mh = np.array(retxt(rep+'Mh_'+str(v0)+'.txt',nbin,0,0))
	refMh, z, xH2r, xHDr, xer, Tbr, vrr = retxt(rep+'ref_'+str(v0)+'.txt',1,0,0)[0]
	plt.figure()
	ax = plt.subplot(111)
	ctf = plt.contourf(X, Y, np.log10(Mh), np.linspace(5.5, 8, 101), cmap=plt.cm.jet)#cm.Blues)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$\log(\tilde{M}_{\mathrm{th}}\ [M_{\odot}])$',size=14)
	cls = [r'$\tilde{M}_{\mathrm{th}}(\mathrm{CDM})$', r'$M_{2}$', r'$10M_{2}$']
	c1 = plt.contour(X, Y, np.log10(Mh), [np.log10(refMh)-1e-2*(v0==0)+1e-2*(v0>40)], colors='k')
	#print(np.min(Mh[Mh!=np.nan]))
	c2 = plt.contour(X, Y, np.log10(Mh), [np.log10(Mup(zvir))], colors='k', linestyles='--')
	c3 = plt.contour(X, Y, np.log10(Mh), [0.99+np.log10(Mup(zvir))], colors='k', linestyles='-.')
	for l in c1.collections:
		l.set_label(cls[0])
	for l in c2.collections:
		l.set_label(cls[1])
	for l in c3.collections:
		l.set_label(cls[2])
	plt.plot([0.3], [8e-20], '*', color='purple', label='Fiducial')
	plt.plot([1e-3], [1e-18], 'o', color='purple', label='Cold')
	plt.plot([10.], [1e-20], '^', color='purple', label='Warm')
	plt.plot(10**llogm, 10**llogs, 'k:', label=r'$\langle\delta T\rangle=-300$ mK')
	if v0==0:
		plt.legend(loc=1)
	ax.tick_params(axis='x', which='minor', bottom=True)
	#lx = X[:,0]
	#plt.plot(lx[lx>10], lx[lx>10]*cons, 'k:')
	#plt.plot([10, 10], [10*cons, np.max(Y)], 'k:')
	plt.xscale('log')
	plt.yscale('log')
	#plt.xlim(np.min(X), np.max(X))
	#plt.ylim(np.min(Y), np.max(Y))
	plt.xticks(np.logspace(-4, 2, 7))
	plt.xlabel(r'$m_{\chi}c^{2}\ [\mathrm{GeV}]$')
	plt.ylabel(r'$\sigma_{1}\ [\mathrm{cm^{2}}]$')
	plt.tight_layout()
	plt.savefig('MthMap_v0_'+str(v0)+'.pdf')
	plt.close()

	d0 = readd(m, zvir, v0, mode = 0, rep = drep)
	d1 = readd(m, zvir, v0, mode = 1, rep = drep)

	lzt = lzt = [10.4]+[20*i+20 for i in range(5)]+[150, 200, 300]#[int(z) for z in np.geomspace(10.0,100.,6)]
	loc = np.log10([TZ(z)/(1e6*YR) for z in lzt])

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	plt.plot(d0['t'], d0['Tb'], label=r'$T_{b}$, CDM', color='b')
	plt.plot(d1['t'], d1['Tb'], '--', label=r'$T_{b}$, BDMS', color='r')
	#plt.plot(d0['t'], d0['Tdm'], '-.', label=r'$T_{\chi}$, CDM')
	plt.plot(d1['t'], d1['Tdm'], '-.', label=r'$T_{\chi}$, BDMS', color='orange')
	plt.plot(d1['t'], (d1['z']+1)*2.726, ':', label=r'$T_{\mathrm{CMB}}$', color='g')
	plt.plot(d1['t'], T_b(d1['z']), 'k-', label=r'$T_{b}(\mathrm{IGM})$, CDM', lw=0.5)
	plt.xlabel(r'$t\ [\mathrm{Myr}]$')
	plt.ylabel(r'$T\ [\mathrm{K}]$')
	if v0==0:
		plt.legend()
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(np.min(d0['t']), np.max(d0['t']))
	plt.ylim(1, 3e3)
	ax2 = ax1.twiny()
	#ax2.set_xscale('log')
	ax2.set_xticks(loc)
	ax2.set_xticklabels([str(int(z)) for z in lzt],size=14)
	ax2.set_xlabel(r'$z$')
	ax2.set_xlim(np.log10(ax1.get_xlim()))
	plt.tight_layout()
	plt.savefig('Example_T_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')
	plt.close()

	if v0==0:
		continue
	vIGM = [vbdm_z(z, v0)/1e5 for z in d0['z']]
	luth = [uthf(PROTON, 0.3*GeV_to_mass, T, Td)/1e5 for T, Td in zip(d1['Tb'], d1['Tdm'])]
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	plt.plot(d0['t'], d0['v']/1e5, label='CDM', color='b')
	plt.plot(d1['t'], d1['v']/1e5, '--', label='BDMS', color='r')
	plt.plot(d1['t'], luth, '-.', label=r'$u_{\mathrm{th}}$', color='orange')
	plt.plot(d0['t'], vIGM, ':', label='IGM, CDM', color='g')
	plt.xlabel(r'$t\ [\mathrm{Myr}]$')
	plt.ylabel(r'$v_{b\chi}\ [\mathrm{km\ s^{-1}}]$')
	plt.legend()
	plt.xlim(np.min(d0['t']), np.max(d0['t']))
	plt.ylim(np.min(vIGM)/10, np.max(vIGM)*10)
	plt.xscale('log')
	plt.yscale('log')
	ax2 = ax1.twiny()
	#ax2.set_xscale('log')
	ax2.set_xticks(loc)
	ax2.set_xticklabels([str(int(z)) for z in lzt],size=14)
	ax2.set_xlabel(r'$z$')
	ax2.set_xlim(np.log10(ax1.get_xlim()))
	plt.tight_layout()
	plt.savefig('Example_vbDM_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')
	plt.close()

v0 = 30
mgas = mmw()*PROTON
nIGM = [rhom(1/(1+z))*0.048/(0.315*mgas) for z in d0['z']]
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(d0['t'], d0['nb'], label='Top-hat model')
	#plt.plot(d1['t'], d1['nb'], '--', label='BDMS')
ax1.plot(d0['t'], nIGM, '-.', label='IGM')
ax1.set_xlabel(r'$t\ [\mathrm{Myr}]$')
ax1.set_ylabel(r'$n\ [\mathrm{cm^{-3}}]$')
ax1.legend()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(np.min(d0['t']), np.max(d0['t']))
ax2 = ax1.twiny()
#lzt = [10]+[20*i+20 for i in range(5)]+[150, 200, 250, 300]#[int(z) for z in np.geomspace(10.0,300.,8)]
#loc = np.log10([TZ(z)/(1e6*YR) for z in lzt])
ax2.set_xlim(np.log10(ax1.get_xlim()))
ax2.set_xticks(loc)
ax2.set_xticklabels([str(int(z)) for z in lzt],size=14)
ax2.set_xlabel(r'$z$')
#ax2.set_xscale('log')
plt.tight_layout()
plt.savefig('Example_n_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')
plt.close()
#"""

lv = np.array(retxt(rep+'vbase.txt',1,0,0)[0])/30 #np.logspace(0, 3, 31)
vd, vu = np.min(lv), np.max(lv)
lz = retxt(rep+'zbase.txt',1,0,0)[0]
nz = len(lz)

lm = np.array(retxt(rep+'Mth_v.txt',nz,0,0))
lxh2 = np.array(retxt(rep+'xh2_v.txt',nz,0,0))
lvr = np.array(retxt(rep+'vr_v.txt',nz,0,0))
lm_ = np.array(retxt(rep+'Mth_v0.txt',nz,0,0))
lxh2_ = np.array(retxt(rep+'xh2_v0.txt',nz,0,0))
lvr_ = np.array(retxt(rep+'vr_v0.txt',nz,0,0))

#"""
y1, y2 = 1e4, 3e9
plt.figure()
a = [plt.plot(lv, lm[i], label=r'$z_{\mathrm{vir}}='+str(int(lz[i]*100)/100)+'$, BDMS', ls = lls[i], color='r', lw = 2) for i in range(nz)]
a = [plt.plot(lv, lm_[i], label=r'$z_{\mathrm{vir}}='+str(int(lz[i]*100)/100)+'$, CDM',ls=lls[i], color='b', lw=1) for i in range(nz)]
plt.plot([1,1],[y1,y2], 'k', lw=0.5)
plt.plot([2,2],[y1,y2], 'k', lw=0.5)
plt.plot([3,3],[y1,y2], 'k', lw=0.5)
plt.legend(loc=2)
plt.xlabel(r'$v_{b\chi,0}\ [\sigma_{\mathrm{rms}}]$')#[\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$\tilde{M}_{\mathrm{th}}\ [M_{\odot}]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(vd, vu)
plt.ylim(y1, y2)
plt.tight_layout()
plt.savefig('Mth_v.pdf')
plt.close()

y1, y2 = 4e-6, 7e-4
plt.figure()
a = [plt.plot(lv, lxh2[i], label=r'$z_{\mathrm{vir}}='+str(int(lz[i]*100)/100)+'$, BDMS', ls = lls[i], color='r', lw = 2) for i in range(nz)]
a = [plt.plot(lv, lxh2_[i], label=r'$z_{\mathrm{vir}}='+str(int(lz[i]*100)/100)+'$, CDM',color='b',ls=lls[i],lw=1) for i in range(nz)]
plt.plot([1,1],[y1,y2], 'k', lw=0.5)
plt.plot([2,2],[y1,y2], 'k', lw=0.5)
plt.plot([3,3],[y1,y2], 'k', lw=0.5)
plt.legend(loc=3)
plt.xlabel(r'$v_{b\chi,0}\ [\sigma_{\mathrm{rms}}]$')#[\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$x_{\mathrm{H_{2}}}$',size=18)
plt.xscale('log')
plt.yscale('log')
plt.xlim(vd, vu)
plt.ylim(y1, y2)
plt.tight_layout()
plt.savefig('XH2_v.pdf')
plt.close()

y1, y2 = 1e-1, 3e3
plt.figure()
a = [plt.plot(lv, lvr[i]/1e5, label=r'$z_{\mathrm{vir}}='+str(int(lz[i]*100)/100)+'$, BDMS', ls = lls[i], color='r', lw=2) for i in range(nz)]
	#a = [plt.plot(lv, [vbdm_z(zr[i], v)/1e5  for v in lv], label=r'$z_{vir}='+str(int(lz[i]*100)/100)+'$, CDM',color='k',ls=lls[i],lw=0.5) for i in range(nz)]
a = [plt.plot(lv, lvr_[i]/1e5, label=r'$z_{\mathrm{vir}}='+str(int(lz[i]*100)/100)+'$, CDM',color='b',ls=lls[i],lw=1) for i in range(nz)]
plt.plot([1,1],[y1,y2], 'k', lw=0.5)
plt.plot([2,2],[y1,y2], 'k', lw=0.5)
plt.plot([3,3],[y1,y2], 'k', lw=0.5)
plt.legend()
plt.xlabel(r'$v_{b\chi,0}\ [\sigma_{\mathrm{rms}}]$')#[\mathrm{km\ s^{-1}}]$')
plt.ylabel(r'$v_{b\chi,V}\ [\mathrm{km\ s^{-1}}]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(vd, vu)
plt.ylim(y1, y2)
plt.tight_layout()
plt.savefig('Vr_v.pdf')
plt.close()
#"""

from matplotlib.ticker import MultipleLocator
nbin = 32
lowb = 1e-7
rep = 'Nhrat/'
lz = [15, 20, 25, 30]
for z in lz:
	d0 = retxt(rep+'nh_ref_z'+str(z)+'.txt',1,0,0)[0]
	X = np.array(retxt(rep+'X_z'+str(z)+'.txt',nbin,0,0))
	Y = np.array(retxt(rep+'Y_z'+str(z)+'.txt',nbin,0,0))
	nh = np.array(retxt(rep+'nh_z'+str(z)+'.txt',nbin,0,0))
	fh = nh/d0
	fh = fh + lowb*(fh<lowb)
	plt.figure()
	ax = plt.subplot(111)
	plt.minorticks_on() 
	ctf = plt.contourf(X, Y, np.log10(fh), np.linspace(-7, 1, 81), cmap=plt.cm.seismic)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$\log(f_{n_{h}})$',size=14)
	ct = plt.contour(X, Y, np.log10(fh), [-3., -2., -1.], colors='k')
	plt.contour(X, Y, np.log10(fh), [0.0-7e-4], linestyles='-', colors='k')
	#plt.clabel(ct, fontsize=12, inline=1)
	plt.plot([0.3], [8e-20], '*', color='purple', label='Fiducial')
	plt.plot([1e-3], [1e-18], 'o', color='purple', label='Cold')
	plt.plot([10.], [1e-20], '^', color='purple', label='Warm')
	plt.plot(10**llogm, 10**llogs, 'k:', label=r'$\langle\delta T\rangle=-300$ mK')
	if z==20:
		plt.legend()
	ml = MultipleLocator(5)
	ax.tick_params(axis='x', which='minor', bottom=True)
	ax.xaxis.set_minor_locator(ml)
	plt.xscale('log')
	plt.yscale('log')
	plt.xticks(np.logspace(-4, 2, 7))
	plt.xlabel(r'$m_{\chi}c^{2}\ [\mathrm{GeV}]$')
	plt.ylabel(r'$\sigma_{1}\ [\mathrm{cm^{2}}]$')
	plt.tight_layout()
	plt.savefig('fhMap_z'+str(z)+'.pdf')
	plt.close()

def Nhrat(rep = 'Nhrat/', mode=0, Mdm = 0.3):
	if mode==1:
		return np.array(retxt(rep+'Nh_ratio.txt',2,0,0))
	lv = retxt(rep+'vbase.txt',1,0,0)[0] #np.logspace(0, 3, 31)
	vd, vu = np.min(lv), np.max(lv)
	lz = np.array(retxt(rep+'zbase.txt',1,0,0)[0])
	nz = len(lz)

	lm = np.array(retxt(rep+'Mth_v.txt',nz,0,0))
	lxh2 = np.array(retxt(rep+'xh2_v.txt',nz,0,0))
	lvr = np.array(retxt(rep+'vr_v.txt',nz,0,0))
	lm_ = np.array(retxt(rep+'Mth_v0.txt',nz,0,0))
	lxh2_ = np.array(retxt(rep+'xh2_v0.txt',nz,0,0))
	lvr_ = np.array(retxt(rep+'vr_v0.txt',nz,0,0))

	lrat = []
	for i in range(nz):
		Nh0 = Nhalo(lz[i], lm_[i], lv, 0, Mdm)
		Nh1 = Nhalo(lz[i], lm[i], lv, 1, Mdm)
		#print('Halo number ratio = {}, at z  = {}'.format(Nh1/Nh0, lz[i]))
		lrat.append(Nh1/Nh0)
	totxt(rep+'Nh_ratio.txt',[lrat, lz],0,0,0)
	lrat = np.array(lrat)
	return [lrat, lz]

lrat, lz = Nhrat('Nhrat/',1)
lrat1, lz1 = Nhrat('Nhrat_1/',1,1e-3)
lrat2, lz2 = Nhrat('Nhrat_2/',1,10.)

ydown, yup = 1e-7, 1e1
x1, x2 = 10, 60
fig = plt.figure()
ax1 = fig.add_subplot(111)
plt.plot(lz[lrat>0], lrat[lrat>0], label='Fiducial')
plt.plot(lz1[lrat1>0], lrat1[lrat1>0], '--', label='Cold')
plt.plot(lz2[lrat2>0], lrat2[lrat2>0], '-.', label='Warm')
plt.xlabel(r'$z_{\mathrm{vir}}$')
plt.ylabel(r'$\bar{n}_{h}^{\mathrm{Pop III}}(\mathrm{BDMS})/\bar{n}_{h}^{\mathrm{Pop III}}(\mathrm{CDM})$')
plt.fill_between([15,20],[ydown,ydown],[yup,yup], facecolor='gray', label='EDGES')
plt.plot([x1, x2], [1, 1], 'k:', lw = 0.5)
#plt.xscale('log')
plt.legend()
plt.yscale('log')
plt.xlim(x1, x2)
plt.ylim(ydown, yup)
ax2 = ax1.twiny()
lzt = np.linspace(x1, x2, 10)
lt = [TZ(z)/(1e6*YR) for z in lzt]
ax2.set_xticks(lzt)
ax2.set_xticklabels([str(int(t*10)/10) for t in lt],size=14)
ax2.set_xlabel(r'$t\ [\mathrm{Myr}]$')
ax2.set_xlim(ax1.get_xlim())
plt.tight_layout()
plt.savefig('Nh_rat.pdf')
#plt.show()
plt.close()

