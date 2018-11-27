from cosmology import *
from bdms import *
from radcool import *
from chemi import *
import os
import multiprocessing as mp
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

def coolt(Tb_old, Tdm_old, v_old, nb, nold, rhob_old, rhodm_old, z, mode, Mdm, sigma, gamma, Om, Ob = 0.048, h = 0.6774, X = 0.76, J_21 = 0, T0 = 2.726):
	xh = 4*X/(1+3*X)
	dTs_dt = [0]*3
	if mode!=0:
		dTs_dt = bdmscool(Tdm_old, Tb_old, v_old, rhob_old, rhodm_old, Mdm, sigma, gamma, X)
	dTb_dt = cool(Tb_old, nb, nold, J_21, z, gamma, X, T0) + dTs_dt[1]
	if dTb_dt == 0:
		return TZ(0)
	else:
		return -Tb_old/dTb_dt
	
def evolve(Mh = 1e6, zvir = 20, z0 = 300, v0 = 30, mode = 0, fac = 1.0, Mdm = 0.3, sigma = 8e-20, num = int(1e3), epsT = 1e-3, epsH = 1e-2, dmax = 18*np.pi**2, gamma = 5/3, X = 0.76, D = 4e-5, Li = 4.6e-10, T0 = 2.726, Om = 0.315, Ob = 0.048, h = 0.6774, dtmin = YR, J_21=0.0, Tmin = 1e-4, vmin = 1e-10, nmax = int(1e6)):
	start = time.time()
	init = initial(z0, v0, mode, Mdm, sigma, T0=T0, Om=Om, Ob=Ob, h=h)
	print(Mdm, sigma, init['Tb'], init['Tdm'], init['vbdm'])
	xh = 4.0*X/(1.0+3.0*X)
	xhe, xd, xli = 1-xh, D, Li
	refa = np.array([xh]*6+[xhe]*3+[xd]*3+[xli]*5)
	t0 = TZ(z0)
	t1 = TZ(zvir)
	tpost = min(fac/H(1/(1+zvir), Om, h), TZ(0)-t1)
	tmax = t1 + tpost
	dt0 = (tmax-t0)/num
	dt0_ = (tmax-t0)/nmax
	#print('Time step: {} yr'.format(dt0/YR))
	rhodm_z = lambda x: rho_z(x, zvir, dmax)*(Om-Ob)/Om
	rhob_z = lambda x: rho_z(x, zvir, dmax)*Ob/Om
	Ns = len(init['X'])
	lz = [z0]
	lt = [t0]
	lTb = [max(init['Tb'], Tmin)]
	lTdm = [max(init['Tdm'], Tmin/1e10)]
	lv = [init['vbdm']]
	lrhob = [rhob_z(z0)]
	lrhodm = [rhodm_z(z0)]
	lX = [[x] for x in init['X']]
	t_cum, count, total = t0, 0, 0
	tag0 = 0
	tag1 = 0
	TV = Tvir(Mh, zvir, dmax)
	VV = Vcir(Mh, zvir, dmax)
	tag2 = 0
	if lTb[0]>TV:
		tag2 = 1
	Tb_V = TV
	pV = []
	para = [Mdm, sigma, gamma, Om, Ob, h, X, J_21, T0]
	tcool = 0.0
	tffV = 1/H(1/(1+zvir), Om, h) #tff(zvir, dmax)
	#yy = np.zeros(Ns, dtype='float')
	tagt = 0
	while t_cum < tmax or z>zvir:
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
				dTs_dt_old = bdmscool(Tdm_old, Tb_old, v_old, rhob_old, rhodm_old, Mdm, sigma, gamma, X)
			dTb_dt_old = cool(Tb_old, nb, nold*nb, J_21, z0, gamma, X, T0) \
						+ dTs_dt_old[1] + dlnrho_dt_old*(gamma-1)*Tb_old \
						+ GammaC(z0, Om, Ob, h = h, X = X, T0 = T0)*(T0*(1+z0)-Tb_old)
			dTdm_dt_old = dTs_dt_old[0] + dlnrho_dt_old*(gamma-1)*Tdm_old
			dv_dt_old = dTs_dt_old[2] + dlnrho_dt_old*(gamma-1)*v_old
		else:
			dt_T = dt0
			if abs(dTb_dt_old*dt_T)>epsT*Tb_old and Tb_old>Tmin*10:
				dt_T = epsT*Tb_old/abs(dTb_dt_old)
				#if dt_T < dt0_ and Tb_old<=Tmin*100:
				#	dt_T = dt0_
		if dt_T + t_cum>t1 and tagt==0:
			dt_T = t1 - t_cum
			tagt = 1
		if dt_T + t_cum>tmax:
			dt_T = tmax - t_cum
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
		uth = (Tb_old*BOL/PROTON+Tdm_old*BOL/(Mdm*GeV_to_mass))**0.5
		if mode!=0 and (Tb_old>Tdm_old or v_old>vmin*uth):
			dTs_dt = bdmscool(Tdm_old, Tb_old, v_old, rhob_old, rhodm_old, Mdm, sigma, gamma, X)
		dTb_dt = cool(Tb_old, nb, nold*nb, J_21, z, gamma, X, T0) + dTs_dt[1] + GammaC(z, Om, Ob, h = h, X = X, T0 = T0)*(T0*(1+z)-Tb_old)
		if tag0==0:
			dTb_dt += dlnrho_dt*(gamma-1)*Tb_old
		dTdm_dt = dTs_dt[0] + dlnrho_dt*(gamma-1)*Tdm_old
		dv_dt = dTs_dt[2] + dlnrho_dt*(gamma-1)*v_old
		Tb_old = max(Tb_old + (dTb_dt + dTb_dt_old)*abund[1]/2.0, Tmin)
		Tdm_old = max(Tdm_old + (dTdm_dt + dTdm_dt_old)*abund[1]/2.0, Tmin/1e10)
		v_old = max(v_old + (dv_dt + dv_dt_old)*abund[1]/2.0, vmin*uth)
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
		if tag2==1:
			if Tb_old<TV:
				tag2 = 0
		#if max(Tb_old, Tdm_old) > TV and tag0==0 and tag2==0:
		#	Tb_old = TV
		#	tag0 = 1
		if t_cum>=t1 and tag1==0 and tag0==0:
			Tb_V = Tb_old
			Tb_old = max(TV, Tb_old)
			Tdm_old = max(TV, Tdm_old)
			#v_old = max(VV, v_old)
			pV = [Tb_old, Tdm_old, v_old, nb, nold*nb, rhob_old, rhodm_old, z]
			tcool = coolt(*pV, mode, *para)
			tag1 = 1
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
	#print(t_cum-t1)
	#print('Time taken: {} s'.format(end-start))
	return d

Mup = lambda z: 2.5e7*((1+z)/10)**-1.5
Mdown = lambda z: 1.54e5*((1+z)/31)**-2.074 #1e6*((1+z)/10)**-2

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

def Mth_z(z1, z2, nzb = 10, m1 = 1e2, m2 = 1e10, nmb = 100, mode = 0, z0 = 300, v0 = 30, Mdm = 0.3, sigma = 8e-20, rat = 1.0, dmax = 18*np.pi**2, Om = 0.315, h = 0.6774, fac = 1e-3):
	m0 = (m1*m2)**0.5
	lz = np.linspace(z1, z2, nzb)
	out = []
	lxh2 = []
	lxhd = []
	for z in lz:
		mmax = Mup(z)*10
		lm = np.logspace(np.log10(m1), np.log10(mmax), nmb)
		d = evolve(m0, z, z0, v0, mode, Mdm = Mdm, sigma = sigma, dmax = dmax, Om = Om, h = h, fac = fac)
		tffV = 1/H(1/(1+z), Om, h) #tff(z, dmax)
		lT = [Tvir(m, z, dmax) for m in lm]
		lv = [Vcir(m, z, dmax) for m in lm]
		print(d['pV'][2])
		pV = d['pV'][2:]
		#pV = d['pV'][3:]
		lt0 = np.array([coolt(T, T, *pV, 1, *d['para'])/tffV for T, v in zip(lT, lv)])
		#lt0 = np.array([coolt(T, T, v, *pV, mode, *d['para'])/tffV for T, v in zip(lT, lv)])
		if np.max(lt0)<=0:
			print('Heating!')
			mth = np.nan
		else:
			lt = lt0[lt0>0]
			lm = lm[lt0>0]
			if np.min(np.log10(lt))>=np.log10(rat):
				mth = np.max(lm)
			elif np.max(np.log10(lt))<=np.log10(rat):
				mth = np.min(lm)
			else:
				rat_m = interp1d(np.log10(lt), np.log10(lm))
				mth = 10**rat_m(np.log10(rat))
		out.append(mth)
		lxh2.append(d['X'][3][-1])
		lxhd.append(d['X'][11][-1])
	return [np.array(out), lz, lxh2, lxhd]

def parasp(v0 = 30., m1 = -4, m2 = 2, s1 = -1, s2 = 4, z = 17, dmax = 200, nbin = 10, fac = 1e-3, rat = 1.0, ncore=4):
	lm = np.logspace(m1, m2, nbin)
	ls = np.logspace(s1, s2, nbin)
	X, Y = np.meshgrid(lm, ls, indexing = 'ij')
	mmax = Mup(17)*10
	lMh = np.zeros(X.shape)
	lXH2 = np.zeros(X.shape)
	lXHD = np.zeros(X.shape)
	np_core = int(nbin/ncore)
	lpr = [[i*np_core, (i+1)*np_core] for i in range(ncore-1)] + [[(ncore-1)*np_core, nbin]]
	print(lpr)
	manager = mp.Manager()
	def sess(pr0, pr1, j):
		out = []
		for i in range(pr0, pr1):
			lm0 = np.logspace(2, np.log10(mmax), 100)
			d = evolve(1e6, z, Mdm = lm[i], sigma = ls[j]*1e-20, v0 = v0, mode = 1, dmax = dmax, fac = fac)
			#print(lm[i], ls[j]*1e-20)
			tffV = 1/H(1/(1+z))
			lT = [Tvir(m, z, dmax) for m in lm0]
			lv = [Vcir(m, z, dmax) for m in lm0]
			print(d['pV'][2])
			pV = d['pV'][2:]
			#pV = d['pV'][3:]
			lt0 = np.array([coolt(T, T, *pV, 1, *d['para'])/tffV for T, v in zip(lT, lv)])
			#lt0 = np.array([coolt(T, T, v, *pV, 1, *d['para'])/tffV for T, v in zip(lT, lv)])
			if np.max(lt0)<=0:
				print('Heating!', np.max(lt0))
				mth = np.nan
			else:
				lt = lt0[lt0>0]
				lm0 = lm0[lt0>0]
				if np.min(np.log10(lt))>=np.log10(rat):
					mth = np.max(lm0)
				elif np.max(np.log10(lt))<=np.log10(rat):
					mth = np.min(lm0)
				else:
					rat_m = interp1d(np.log10(lt), np.log10(lm0))
					mth = 10**rat_m(np.log10(rat))
			out.append([mth, d['X'][3][-1], d['X'][11][-1]])
		output.put((pr0, np.array(out).T))
	for i in range(nbin):
		output = manager.Queue()
		pros = [mp.Process(target=sess, args=(lpr[k][0], lpr[k][1], i)) for k in range(ncore)]
		for p in pros:
			p.start()
		for p in pros:
			p.join()
		out = [output.get() for p in pros]
		out.sort()
		lMh[:,i] = np.hstack([x[1][0] for x in out])
		lXH2[:,i] = np.hstack([x[1][1] for x in out])
		lXHD[:,i] = np.hstack([x[1][2] for x in out])
		#for j in range(nbin):
		#	sol = T21_pred(v0, lm[i], ls[j]*1e-20, xa0)
		#	lT[i,j] = -sol[0]
		#	lTb[i,j] = sol[1]
	return X, Y*1e-20, lMh, lXH2, lXHD

if __name__=="__main__":
	tag = 0
	v0 = 30
	nbin = 32
	ncore = 4
	dmax = 18*np.pi**2 * 1
	rat = 1.
	fac = 1e-3
	rep = '1-sigma_ns/'
	if not os.path.exists(rep):
		os.makedirs(rep)

	#refMh = Mth_z(16, 17, 2, v0 = v0, dmax = dmax, Mdm = 0.3, sigma = 8e-20, mode = 0, rat = rat)[0][1]
	#refd = stored(evolve(1e6, 17, v0 = v0, fac = fac, mode = 0, dmax = dmax), 1e6, 17, 0)
	#xH2r = refd['X'][3][-1]
	#xHDr = refd['X'][11][-1]
	#totxt(rep+'ref.txt', [[refMh, xH2r, xHDr]], 0,0,0)

	#"""
	if tag==0:
		refMh = Mth_z(16, 17, 2, v0 = v0, dmax = dmax, Mdm = 0.3, sigma = 8e-20, mode = 0, rat = rat)[0][1]
		refd = stored(evolve(1e6, 17, v0 = v0, fac = fac, mode = 0, dmax = dmax), 1e6, 17, 0)
		xH2r = refd['X'][3][-1]
		xHDr = refd['X'][11][-1]
		totxt(rep+'ref.txt', [[refMh, xH2r, xHDr]], 0,0,0)
		X, Y, Mh, XH2, XHD = parasp(v0, m1 = -4, m2 = 2, s1 = -1, s2 = 4, nbin = nbin, ncore = ncore, dmax = dmax, fac = fac, rat = rat)
		totxt(rep+'X_'+str(v0)+'.txt',X,0,0,0)
		totxt(rep+'Y_'+str(v0)+'.txt',Y,0,0,0)
		totxt(rep+'Mh_'+str(v0)+'.txt',Mh,0,0,0)
		totxt(rep+'XH2_'+str(v0)+'.txt',XH2,0,0,0)
		totxt(rep+'XHD_'+str(v0)+'.txt',XHD,0,0,0)
	else:
		X = np.array(retxt(rep+'X_'+str(v0)+'.txt',nbin,0,0))
		Y = np.array(retxt(rep+'Y_'+str(v0)+'.txt',nbin,0,0))
		Mh = np.array(retxt(rep+'Mh_'+str(v0)+'.txt',nbin,0,0))
		XH2 = np.array(retxt(rep+'XH2_'+str(v0)+'.txt',nbin,0,0))
		XHD = np.array(retxt(rep+'XHD_'+str(v0)+'.txt',nbin,0,0))

	refMh, xH2r, xHDr = retxt(rep+'ref.txt',1,0,0)[0]
	print('Reference mass thresold: {} 10^6 Msun'.format(refMh/1e6))
	print('Reference H2 abundance: {} * 10^-4'.format(xH2r*1e4))
	print('Reference HD abundance: {} * 10^-3'.format(xHDr*1e3))
	Mbd = Mup(17)*2
	Mh = Mh*(Mh<Mbd) + Mbd*(Mh>=Mbd)
	plt.figure()
	ctf = plt.contourf(X, Y, np.log10(Mh/1e6), 2*nbin, cmap=plt.cm.Blues)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$\log(M_{\mathrm{th}}\ [10^{6}\ M_{\odot}])$',size=12)
	plt.contour(X, Y, np.log10(Mh/1e6), [np.log10(refMh/1e6)], colors='k')
	print(np.min(Mh[Mh!=np.nan]/1e6))
	plt.contour(X, Y, np.log10(Mh/1e6), [np.log10(Mup(17)/1e6)], colors='k', linestyles='--')
	plt.plot([0.3], [8e-20], '*', color='purple')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$m_{\mathrm{DM}}c^{2}\ [\mathrm{GeV}]$')
	plt.ylabel(r'$\sigma_{1}\ [\mathrm{cm^{2}}]$')
	plt.tight_layout()
	plt.savefig(rep+'MthMap_v0_'+str(v0)+'.pdf')

	plt.figure()
	ctf = plt.contourf(X, Y, XH2*1e4, 2*nbin, cmap=plt.cm.Blues)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$x_{\mathrm{H_{2}}}\cdot 10^{4}$',size=12)
	plt.contour(X, Y, XH2*1e4, [xH2r*1e4], colors='k')
	plt.plot([0.3], [8e-20], '*', color='purple')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$m_{\mathrm{DM}}c^{2}\ [\mathrm{GeV}]$')
	plt.ylabel(r'$\sigma_{1}\ [\mathrm{cm^{2}}]$')
	plt.tight_layout()
	plt.savefig(rep+'XH2Map_v0_'+str(v0)+'.pdf')

	plt.figure()
	ctf = plt.contourf(X, Y, XHD*1e3, 2*nbin, cmap=plt.cm.Blues)
	for c in ctf.collections:
		c.set_edgecolor('face')
	cb = plt.colorbar()
	cb.set_label(r'$x_{\mathrm{HD}}\cdot 10^{3}$',size=12)
	plt.contour(X, Y, XHD*1e3, [xHDr*1e3], colors='k')
	plt.plot([0.3], [8e-20], '*', color='purple')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$m_{\mathrm{DM}}c^{2}\ [\mathrm{GeV}]$')
	plt.ylabel(r'$\sigma_{1}\ [\mathrm{cm^{2}}]$')
	plt.tight_layout()
	plt.savefig(rep+'XHDMap_v0_'+str(v0)+'.pdf')
	#"""

	#"""

	#lm, lz, lxh2, lxhd = Mth_z(10, 100, 46, mode = 0, rat = rat, dmax = dmax, fac = fac)
	#totxt(rep+'Mthz_CDM.txt',[lz, lm, lxh2, lxhd],0,0,0)

	tag = 0
	#v0 = 0.1
	#rat = 10.
	if tag==0:
		lm_, lz_, lxh2_, lxhd_ = Mth_z(10, 100, 46, mode = 1, v0 = v0, rat = rat, dmax = dmax, fac = fac)
		lm, lz, lxh2, lxhd = Mth_z(10, 100, 46, mode = 0, rat = rat, dmax = dmax, fac = fac)
		totxt(rep+'Mthz_CDM.txt',[lz, lm, lxh2, lxhd],0,0,0)
		totxt(rep+'Mthz_BDMS_'+str(v0)+'.txt',[lz_, lm_, lxh2_, lxhd_],0,0,0)
	lz_, lm_, lxh2_, lxhd_ = np.array(retxt(rep+'Mthz_BDMS_'+str(v0)+'.txt',4,0,0))
	lz, lm, lxh2, lxhd = np.array(retxt(rep+'Mthz_CDM.txt',4,0,0))
	plt.figure()
	plt.plot(lz, lm, label='CDM')
	plt.plot(lz_, lm_, '--', label='BDMS')
	plt.plot(lz, Mdown(lz), 'k-.', label='Trenti & Stiavelli (2009)')
	plt.legend()
	plt.xlabel(r'$z_{\mathrm{vir}}$')
	plt.ylabel(r'$M_{\mathrm{th}}\ [M_{\odot}]$')
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'Mth_z_'+str(v0)+'.pdf')
	plt.figure()
	plt.plot(lz, lxh2, label='CDM')
	plt.plot(lz_, lxh2_, '--', label='BDMS')
	plt.legend()
	plt.xlabel(r'$z_{\mathrm{vir}}$')
	plt.ylabel(r'$x_{\mathrm{H_{2}}}$')
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'XH2_z_'+str(v0)+'.pdf')
	plt.figure()
	plt.plot(lz, lxhd, label='CDM')
	plt.plot(lz_, lxhd_, '--', label='BDMS')
	plt.legend()
	plt.xlabel(r'$z_{\mathrm{vir}}$')
	plt.ylabel(r'$x_{\mathrm{HD}}$')
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep+'XHD_z_'+str(v0)+'.pdf')
	#plt.show()
	#"""
	
	#"""
	lls = ['-', '--', '-.', ':']*2
	tag = 1
	#rat = 10.
	nz = 3
	z1, z2 = 20, 100

	#d0 = Mth_z(z1, z2, nz, mode = 0, rat = rat, dmax = dmax, fac = fac)
	#totxt(rep+'ref_z.txt', d0, 0, 0)

	lv = np.logspace(-1, 3, 50)
	vd, vu = np.min(lv), np.max(lv)
	if tag==0:
		d0 = Mth_z(z1, z2, nz, mode = 0, rat = rat, dmax = dmax, fac = fac)
		totxt(rep+'ref_z.txt', d0, 0, 0)
		out = []
		for v in lv:
			d = Mth_z(z1, z2, nz, mode = 1, v0 = v, rat = rat, dmax = dmax, fac = fac)
			out.append(d)
		lz = d[1]
		lm = [[x[0][i] for x in out] for i in range(len(d[0]))]
		lxh2 = [[x[2][i] for x in out] for i in range(len(d[0]))]
		lxhd = [[x[3][i] for x in out] for i in range(len(d[0]))]
		totxt(rep+'Mth_v.txt',lm,0,0,0)
		totxt(rep+'xh2_v.txt',lxh2,0,0,0)
		totxt(rep+'xhd_v.txt',lxhd,0,0,0)
		totxt(rep+'zbase.txt',[lz],0,0,0)
		totxt(rep+'vbase.txt',[lv],0,0,0)
	d0 = retxt(rep+'ref_z.txt',4,0,0)
	mr, zr, xh2r, xhdr = d0
	lv = np.array(retxt(rep+'vbase.txt',1,0,0)[0])
	lz = retxt(rep+'zbase.txt',1,0,0)[0]
	nz = len(lz)
	lm = np.array(retxt(rep+'Mth_v.txt',nz,0,0))
	lxh2 = np.array(retxt(rep+'xh2_v.txt',nz,0,0))
	lxhd = np.array(retxt(rep+'xhd_v.txt',nz,0,0))
	plt.figure()
	a = [plt.plot(lv, lm[i], label=r'$z='+str(int(lz[i]*100)/100)+'$', ls = lls[i]) for i in range(nz)]
	a = [plt.plot(lv, np.ones(len(lv))*mr[i], label=r'$z='+str(int(lz[i]*100)/100)+'$, CDM',color='k',ls=lls[i],lw=0.5) for i in range(nz)]
	plt.legend()
	plt.xlabel(r'$v_{\mathrm{bDM},0}\ [\mathrm{km\ s^{-1}}]$')
	plt.ylabel(r'$M_{\mathrm{th}}\ [M_{\odot}]$')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(vd, vu)
	plt.tight_layout()
	plt.savefig(rep+'Mth_v.pdf')
	plt.figure()
	a = [plt.plot(lv, lxh2[i], label=r'$z='+str(int(lz[i]*100)/100)+'$', ls = lls[i]) for i in range(nz)]
	a = [plt.plot(lv, np.ones(len(lv))*xh2r[i], label=r'$z='+str(int(lz[i]*100)/100)+'$, CDM',color='k',ls=lls[i],lw=0.5) for i in range(nz)]
	plt.legend()
	plt.xlabel(r'$v_{\mathrm{bDM},0}\ [\mathrm{km\ s^{-1}}]$')
	plt.ylabel(r'$x_{\mathrm{H_{2}}}$')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(vd, vu)
	plt.tight_layout()
	plt.savefig(rep+'XH2_v.pdf')
	plt.figure()
	a = [plt.plot(lv, lxhd[i], label=r'$z='+str(int(lz[i]*100)/100)+'$', ls = lls[i]) for i in range(nz)]
	a = [plt.plot(lv, np.ones(len(lv))*xhdr[i], label=r'$z='+str(int(lz[i]*100)/100)+'$, CDM',color='k',ls=lls[i],lw=0.5) for i in range(nz)]
	plt.legend()
	plt.xlabel(r'$v_{\mathrm{bDM},0}\ [\mathrm{km\ s^{-1}}]$')
	plt.ylabel(r'$x_{\mathrm{HD}}$')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(vd, vu)
	plt.tight_layout()
	plt.savefig(rep+'XHD_v.pdf')
	#"""

	"""
	tag = 0
	m = 1e6
	zvir = 10
	v0 = 30
	rep0 = 'example/'
	dmax = 18 * np.pi**2 * 1
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
	plt.plot(d0['t'], d0['Tb'], label=r'$T_{\mathrm{b}}$, CDM')
	plt.plot(d1['t'], d1['Tb'], '--', label=r'$T_{\mathrm{b}}$, BDMS')
	#plt.plot(d0['t'], d0['Tdm'], '-.', label=r'$T_{\mathrm{DM}}$, CDM')
	plt.plot(d1['t'], d1['Tdm'], '-.', label=r'$T_{\mathrm{DM}}$, BDMS')
	plt.xlabel(r'$t\ [\mathrm{yr}]$')
	plt.ylabel(r'$T\ [\mathrm{K}]$')
	plt.legend()
	plt.xscale('log')
	plt.yscale('log')
	plt.tight_layout()
	plt.savefig(rep0+'Example_T_t_m6_'+str(m/1e6)+'_z_'+str(zvir)+'_v0_'+str(v0)+'.pdf')

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


