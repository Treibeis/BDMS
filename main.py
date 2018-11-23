from cosmology import *
from bdms import *
from radcool import *
from chemi import *
import time

x0_default = [1., 5e-4, 1e-18, 1e-11, 5e-16] + \
			 [5e-4, 1.0, 1e-19, 0.0] + \
			 [1.0, 5e-4, 2.5e-11] + \
			 [1.0, 0, 0, 0, 0]

def initial(z0 = 300, v0 =30, mode = 0, Mdm = 0.3, sigma = 8e-20, x0 = x0_default, z00 = 1000, Om = 0.315, Ob = 0.048, h = 0.6774, vmin = 1e-10, T0 = 2.726):
	#0: XH, 1: XH+, 2: XH-, 3: XH2, 4: XH2+,
	#5: Xe, 6: XHe, 7: XHe+, 8: XHe++, 
	#9: XD, 10: XD+, 11: XHD,
	#12: Li, 13: Li+, 14: Li-, 15: LiH+, 16: LiH
	d = {}
	if mode!=0:
		d0 = thermalH(z00, z0, v0, Mdm, sigma, Om, Ob, T0=T0, h=h)
		d['Tb'] = d0['Tb'][-1]
		d['Tdm'] = d0['Tdm'][-1]
		d['vbdm'] = max(d0['v'][-1], vmin)
	else:
		d['Tb'] = T_b(z0, T0=T0)
		d['Tdm'] = T_dm(z0, Mdm, T0=T0)
		d['vbdm'] = vbdm_z(z0, v0)
	d['X']  = x0
	return d

def M_T(T, z, delta, Om = 0.315):
	return (T/(delta*Om/(200*0.315))**(1/3)*10/(1+z)/190356.40337133306)**(3/2)*1e10

def coolt(Tb_old, Tdm_old, nb, nold, v_old, rhob_old, rhodm_old, z, mode, Mdm, sigma, gamma, Om, Ob, h, X, J_21, T0):
	xh = 4*X/(1+3*X)
	dTs_dt = [0]*3
	if mode!=0:
		dTs_dt = bdmscool(Tdm_old, Tb_old, v_old, z, rhob_old, rhodm_old, Mdm, sigma, gamma, X, Om, Ob, h)
	dTb_dt = cool(Tb_old, nb, nold, J_21, z, gamma, X, T0) + dTs_dt[1]
	return -Tb_old/dTb_dt
	
def evolve(Mh = 1e6, zvir = 20, z0 = 300, v0 = 30, mode = 0, fac = 0.1, Mdm = 0.3, sigma = 8e-20, num = int(1e3), epsT = 1e-3, epsH = 1e-2, dmax = 18*np.pi**2, gamma = 5/3, X = 0.76, D = 4e-5, Li = 4.6e-10, T0 = 2.726, Om = 0.315, Ob = 0.048, h = 0.6774, dtmin = 1e3*YR, J_21=0.0, Tmin = 1e-10, vmin = 1e-10):
	start = time.time()
	init = initial(z0, v0, mode, Mdm, sigma, T0=T0, Om=Om, Ob=Ob, h=h)
	xh = 4.0*X/(1.0+3.0*X)
	xhe, xd, xli = 1-xh, D, Li
	refa = np.array([xh]*6+[xhe]*3+[xd]*3+[xli]*5)
	t0 = TZ(z0)
	t1 = TZ(zvir)
	tpost = min(fac/H(1/(1+zvir), Om, h), TZ(0)-t1)
	tmax = t1 + tpost
	dt0 = (tmax-t0)/num
	print('Time step: {} yr'.format(dt0/YR))
	rhodm_z = lambda x: rho_z(x, zvir, dmax)*(Om-Ob)/Om
	rhob_z = lambda x: rho_z(x, zvir, dmax)*Ob/Om
	Ns = len(init['X'])
	lz = [z0]
	lt = [t0]
	lTb = [init['Tb']]
	lTdm = [init['Tdm']]
	lv = [init['vbdm']]
	lrhob = [rhob_z(z0)]
	lrhodm = [rhodm_z(z0)]
	lX = [[x] for x in init['X']]
	t_cum, count, total = t0, 0, 0
	tag0 = 0
	tag1 = 0
	TV = Tvir(Mh, zvir, dmax)
	tag2 = 0
	if lTb[0]>TV:
		tag2 = 1
	Tb_V = TV
	pV = []
	para = [Mdm, sigma, gamma, Om, Ob, h, X, J_21, T0]
	tcool = 0.0
	tffV = 1/H(1/(1+zvir), Om, h) #tff(zvir, dmax)
	#yy = np.zeros(Ns, dtype='float')
	while t_cum < tmax:
		if count==0:
			z = z0
			dt_T = dtmin
			yy = np.array([x[count] for x in lX])
			mgas = mmw(yy[5], yy[7], X)*PROTON
			nb = lrhob[0]/mgas
			nold = yy * refa
			Tb_old = lTb[count]
			Tdm_old = lTdm[count]
			v_old = lv[count]
			rhob_old, rhodm_old = lrhob[count], lrhodm[count]
			dlnrho_dt_old = Dlnrho(t0, t0+dtmin/2.0, zvir, dmax)
			dTs_dt_old = [0]*3
			dTs_dt = [0]*3
			if mode!=0:
				dTs_dt_old = bdmscool(Tdm_old, Tb_old, v_old, z0, rhob_old, rhodm_old, Mdm, sigma, gamma, X, Om, Ob, h, T0 = T0)
			dTb_dt_old = cool(Tb_old, nb, nold*nb, J_21, z0, gamma, X, T0) + dTs_dt_old[1] + dlnrho_dt_old*(gamma-1)*Tb_old
			dTdm_dt_old = dTs_dt_old[0] + dlnrho_dt_old*(gamma-1)*Tdm_old
			dv_dt_old = dTs_dt_old[2] + dlnrho_dt_old*(gamma-1)*v_old
		else:
			dt_T = dt0
			if abs(dTb_dt_old*dt_T)>epsT*Tb_old:
				dt_T = max(epsT*Tb_old/abs(dTb_dt_old), dtmin)
		if dt_T + t_cum>tmax:
			dt_T = tmax - t_cum
			t_cum = tmax
		if tag2==1:
			if Tb_old<TV:
				tag2 = 0
		#if max(Tb_old, Tdm_old) > TV and tag0==0 and tag2==0:
		#	Tb_old = TV
		#	tag0 = 1
		if z<zvir and tag1==0 and tag0==0:
			Tb_V = Tb_old
			Tb_old = max(TV, Tb_old)
			Tdm_old = max(TV, Tdm_old)
			pV = [Tb_old, Tdm_old, nb, nold*nb, v_old, rhob_old, rhodm_old, z]
			tcool = coolt(*pV, mode, *para)
			tag1 = 1
		if count==0:
			Cr0, Ds0 = np.zeros(Ns,dtype='float'), np.zeros(Ns,dtype='float')
			abund0 = chemistry1(Tb_old, nold*nb, dt_T, epsH, J_21, Ns, xh*nb, xhe*nb, xd*nb, xli*nb, Cr0, Ds0)
			Cr0, Ds0 = abund0[5], abund0[6]
		else:
			Cr0, Ds0 = abund[5], abund[6]
		nold = yy * refa
		abund = chemistry1(Tb_old, nold*nb, dt_T, epsH, J_21, Ns, xh*nb, xhe*nb, xd*nb, xli*nb, Cr0, Ds0)
		nold = abund[0]/nb
		for x in range(Ns):
			if refa[x]!=0:
				yy[x] = nold[x]/refa[x]
			else:
				yy[x] = 0.0
		mgas = mmw(yy[5], yy[7], X)*PROTON
		#if count<10:
		#	print(nold, Tb_old)
		t_cum += abund[1]
		z = z_t(t_cum)
		dlnrho_dt = Dlnrho(t_cum, t_cum + abund[1]/2.0, zvir, dmax)
		if mode!=0:
			dTs_dt = bdmscool(Tdm_old, Tb_old, v_old, z, rhob_old, rhodm_old, Mdm, sigma, gamma, X, Om, Ob, h, T0 = T0)
		dTb_dt = cool(Tb_old, nb, nold*nb, J_21, z, gamma, X, T0) + dTs_dt[1]
		if tag0==0:
			dTb_dt += dlnrho_dt*(gamma-1)*Tb_old
		dTdm_dt = dTs_dt[0] + dlnrho_dt*(gamma-1)*Tdm_old
		dv_dt = dTs_dt[2] + dlnrho_dt*(gamma-1)*v_old
		Tb_old = max(Tb_old + (dTb_dt + dTb_dt_old)*abund[1]/2.0, Tmin)
		Tdm_old = max(Tdm_old + (dTdm_dt + dTdm_dt_old)*abund[1]/2.0, Tmin)
		v_old = max(v_old + (dv_dt + dv_dt_old)*abund[1]/2.0, vmin)
		dTb_dt_old = dTb_dt
		dTdm_dt_old = dTdm_dt
		dv_dt_old = dv_dt
		if tag0==0:
			rhob_old = rhob_z(z)
		rhodm_old = rhodm_z(z)
		nb = rhob_old/mgas
		#print(z, nb, mgas)
		total += abund[4]
		count += 1
		if (count%10==0)or(t_cum>=tmax):
			lt.append(t_cum)#[count] = t_cum
			lTb.append(Tb_old)#[count] = Told
			lTdm.append(Tdm_old)#[count] = nb
			lv.append(v_old)
			lrhob.append(rhob_old)
			lrhodm.append(rhodm_old)
			lz.append(z)
			for x in range(Ns):
				if refa[x]!=0:
					lX[x].append(nold[x]/refa[x])#[count] = newX[x]
				else:
					lX[x].append(0.0)
	d = {}
	d['t'] = np.array(lt)/YR/1e6
	d['z'] = np.array(lz)
	d['Tb'] = np.array(lTb)
	d['Tdm'] = np.array(lTdm)
	d['v'] = np.array(lv)
	d['rho'] = np.array(lrhob) + np.array(lrhodm)
	d['nb'] = np.array(lrhob)/mgas
	d['X'] = np.array(lX)
	d['rat'] = Tb_old/TV
	d['rat0'] = tpost/(t1 + tpost)
	d['s'] = int(tpost/tmax > Tb_old/TV)
	d['Tvir'] = TV
	d['TbV'] = Tb_V
	d['rat1'] = Tb_V/TV
	d['rat2'] = tcool/tffV
	d['m'] = M_T(Tb_V/d['rat0'], zvir, dmax)
	d['pV'] = pV
	d['para'] = para
	end = time.time()
	print('Time taken: {} s'.format(end-start))
	return d

def stored(d, Mh = 1e6, zvir = 20, v0 = 30, mode = 0, Mdm = 0.3, sigma = 8e-20, rep = 'data/'):
	out0 = [d['t'], d['z'], d['Tb'], d['Tdm'], d['v'], d['rho'], d['nb']]
	out1 = [[d['rat'], d['rat0'], d['rat1'], d['rat2'], d['s'], d['Tvir'], d['TbV'], d['m']]]
	base = 'M'+str(int(Mh/1e6 * 100)/100)+'_z'+str(zvir)+'_v'+str(int(v0*100)/100)
	if mode!=0:
		base = base + '_Mdm'+str(Mdm)+'_sigma'+str(sigma)
	totxt(rep+'dataD_'+base+'.txt',out0,0,0,0)
	totxt(rep+'dataX_'+base+'.txt',d['X'],0,0,0)
	totxt(rep+'dataP_'+base+'.txt',out1,0,0,0)
	return d

def readd(Mh = 1e6, zvir = 20, v0 = 30, mode = 0, Mdm = 0.3, sigma = 8e-20, dmax = 18*np.pi**2, rep = 'data/'):
	base = 'M'+str(int(Mh/1e6 * 100)/100)+'_z'+str(zvir)+'_v'+str(int(v0*100)/100)
	if mode!=0:
		base = base + '_Mdm'+str(Mdm)+'_sigma'+str(sigma)
	rd0 = np.array(retxt(rep+'dataD_'+base+'.txt',7,0,0))
	rd1 = np.array(retxt(rep+'dataP_'+base+'.txt',1,0,0)[0])
	d = {}
	d['X'] = np.array(retxt(rep+'dataX_'+base+'.txt',17,0,0))
	d['t'], d['z'], d['Tb'], d['Tbm'], d['v'], d['rho'], d['nb'] = \
		rd0[0], rd0[1], rd0[2], rd0[3], rd0[4], rd0[5], rd0[6]
	d['rat'], d['rat0'], d['rat1'], d['rat2'], d['s'], d['Tvir'], d['TbV'] = \
		rd1[0], rd1[1], rd1[2], rd1[3], rd1[4], rd1[5], rd1[6]
	d['m'] = M_T(d['TbV']/d['rat0'], zvir, dmax)
	return d

def Mth_z(z1, z2, nzb = 10, m1 = 1e4, m2 = 1e8, nmb = 100, mode = 0, z0 = 300, v0 = 30, Mdm = 0.3, sigma = 8e-20, rat = 1.0, dmax = 18*np.pi**2, Om = 0.315, h = 0.6774):
	m0 = (m1*m2)**0.5
	lm = np.logspace(np.log10(m1), np.log10(m2), nmb)
	lz = np.linspace(z1, z2, nzb)
	out = []
	for z in lz:
		d = evolve(m0, z, z0, v0, mode, Mdm = Mdm, sigma = sigma, dmax = dmax, Om = Om, h = h)
		tffV = 1/H(1/(1+z), Om, h) #tff(z, dmax)
		lT = [Tvir(m, z, dmax) for m in lm]
		pV = d['pV'][2:]
		lt = [max(coolt(T, T, *pV, mode, *d['para'])/tffV, 1e-4) for T in lT]
		#print(lt)
		rat_m = interp1d(np.log10(lt), np.log10(lm))
		mth = 10**rat_m(np.log10(rat))
		#print(mth)
		out.append(mth)
	return np.array(out), lz

if __name__=="__main__":
	tag = 0
	m = 1e6
	zvir = 70
	v0 = 30
	rep0 = 'example/'
	dmax = 200 #1e3
	if tag==0:
		d0 = stored(evolve(m, zvir, mode = 0, dmax = dmax, v0 = v0), m, zvir, mode = 0)
		d1 = stored(evolve(m, zvir, mode = 1, dmax = dmax, v0 = v0), m, zvir, mode = 1)
	else:
		d0 = readd(m, zvir, mode = 0)
		d1 = readd(m, zvir, mode = 1)
	
	mgas = mmw()*PROTON
	nIGM = [rhom(1/(1+z))*0.048/(0.315*mgas) for z in d0['z']]
	plt.figure()
	plt.plot(d0['t'], d0['nb'], label='CDM')
	plt.plot(d1['t'], d1['nb'], '--', label='BDMS')
	plt.plot(d0['t'], nIGM, '-.', label='IGM')
	plt.xlabel(r'$t\ [\mathrm{yr}]$')
	plt.ylabel(r'$n\ [\mathrm{cm^{-3}}]$')
	plt.legend()
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep0+'Example_n_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')

	vIGM = [vbdm_z(z, v0)/1e5 for z in d0['z']]
	plt.figure()
	plt.plot(d0['t'], d0['v']/1e5, label='CDM')
	plt.plot(d1['t'], d1['v']/1e5, '--', label='BDMS')
	plt.plot(d0['t'], vIGM, '-.', label='IGM')
	plt.xlabel(r'$t\ [\mathrm{yr}]$')
	plt.ylabel(r'$v_{\mathrm{bDM}}\ [\mathrm{km\ s^{-1}}]$')
	plt.legend()
	plt.ylim(np.min(vIGM)/10, np.max(vIGM)*10)
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep0+'Example_vdDM_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')

	plt.figure()
	plt.plot(d0['t'], d0['Tb'], label='CDM')
	plt.plot(d1['t'], d1['Tb'], '--', label='BDMS')
	plt.xlabel(r'$t\ [\mathrm{yr}]$')
	plt.ylabel(r'$T_{\mathrm{b}}\ [\mathrm{K}]$')
	plt.legend()
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep0+'Example_Tb_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')

	plt.figure()
	plt.plot(d0['t'], d0['X'][3], label='CDM')
	plt.plot(d1['t'], d1['X'][3], '--', label='BDMS')
	plt.xlabel(r'$t\ [\mathrm{yr}]$')
	plt.ylabel(r'$x_{\mathrm{H_{2}}}$')
	plt.legend()
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim(1e-5, np.max(d0['X'][3])*1.5)
	plt.tight_layout()
	plt.savefig(rep0+'Example_xH2_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')

	plt.figure()
	plt.plot(d0['t'], d0['X'][5], label='CDM')
	plt.plot(d1['t'], d1['X'][5], '--', label='BDMS')
	plt.xlabel(r'$t\ [\mathrm{yr}]$')
	plt.ylabel(r'$x_{\mathrm{e}}$')
	plt.legend()
	plt.xscale('log')
	plt.yscale('log')
	#plt.ylim(1e-5, 1e-3)
	plt.tight_layout()
	plt.savefig(rep0+'Example_xe_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'.pdf')

	plt.figure()
	plt.plot(d0['t'], d0['X'][0], label='0')
	plt.plot(d0['t'], d0['X'][3], label='3')
	plt.plot(d0['t'], d0['X'][5], label='5')
	plt.plot(d0['t'], d0['X'][9], label='9')
	plt.plot(d0['t'], d0['X'][10], label='10')
	plt.plot(d0['t'], d0['X'][11], label='11')
	#plt.plot(d1['t'], d1['X'][5], '--', label='BDMS')
	plt.xlabel(r'$t\ [\mathrm{yr}]$')
	plt.ylabel(r'$x$')
	plt.legend()
	plt.xscale('log')
	plt.yscale('log')
	#plt.ylim(1e-5, 1e-3)
	plt.tight_layout()
	plt.savefig(rep0+'Example_X_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'.pdf')
	#plt.show()

	"""
	#mode = 0
	Mdm, sigma = 0.3, 1e-19
	lm_, lz_ = Mth_z(10, 100, 19, mode = 1)
	lm, lz = Mth_z(10, 100, 19, mode = 0)
	plt.figure()
	plt.plot(lz, lm, label='CDM')
	plt.plot(lz_, lm_, '--', label='BDMS')
	plt.legend()
	plt.xlabel(r'$z$')
	plt.ylabel(r'$M_{\mathrm{th}}\ [M_{\odot}]$')
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig('Mth_z.pdf')
	totxt('Mthz_CDM.txt',[lz, lm],0,0,0)
	totxt('Mthz_BDMS.txt',[lz_, lm_],0,0,0)
	#plt.show()
	"""

	"""
	tag = 0
	zvir = 50
	lm = 10**np.linspace(4,8,10)
	if tag==0:
		ld0 = [stored(evolve(m, zvir), m, zvir) for m in lm]
		#a = [stored(d, m, zvir) for d, m in zip(ld0, lm)]
	else:
		ld0 = [readd(m, zvir) for m in lm]
	lrat = np.array([d['rat'] for d in ld0])
	lrat0 = np.array([d['rat0'] for d in ld0])
	lrat1 = np.array([d['rat1'] for d in ld0])
	lrat2 = np.array([d['rat2'] for d in ld0])
	
	plt.figure()
	plt.plot(lm, lrat/lrat0, label=r'$(T_{f}/T_{V})/(\Delta t/t_{f})$')
	plt.plot(lm, lrat1/lrat0, '--', label=r'$(T(z_{V})/T_{V})/(\Delta t/t_{f})$')
	plt.plot(lm, lrat2, '-.', label=r'$t_{\mathrm{cool}}/t_{\mathrm{ff}}$')
	plt.plot(lm, np.ones(len(lm)), 'k:')
	plt.legend()
	plt.xlim(np.min(lm), np.max(lm))
	plt.ylim(1e-2, 1e2)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$M\ [M_{\odot}]$')
	plt.ylabel('Ratio')
	plt.tight_layout()
	plt.savefig('Ratio_z'+str(zvir)+'.pdf')
	#plt.show()

	mth = [d['m'] for d in ld0]
	print(mth)
	"""


