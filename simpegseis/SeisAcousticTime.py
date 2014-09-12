from SimPEG import *
from SimPEG.Utils import sdiag, mkvc, sdInv
import matplotlib.pyplot as plt
from time import clock

class AcousticTx(Survey.BaseTx):


	def __init__(self, loc, time, rxList, **kwargs):

		self.dt = time[1]-time[0]
		self.time = time
		self.loc = loc
		self.rxList = rxList
		self.kwargs = kwargs


	def RickerWavelet(self):

		"""
			Generating Ricker Wavelet

			.. math ::

		"""
		tlag = self.kwargs['tlag']
		fmain = self.kwargs['fmain']
		time = self.time
		self.wave = (1-2*np.pi**2*fmain**2*(time-tlag)**2)*np.exp(-np.pi**2*fmain**2*(time-tlag)**2)
		return self.wave

	def Wave(self, tInd):

		"""
			Generating Ricker Wavelet

			.. math ::

		"""
		tlag = self.kwargs['tlag']
		fmain = self.kwargs['fmain']
		time = self.time[tInd]
		self.wave = (1-2*np.pi**2*fmain**2*(time-tlag)**2)*np.exp(-np.pi**2*fmain**2*(time-tlag)**2)
		return self.wave

	def getq(self, mesh):

		txind = Utils.closestPoints(mesh, self.loc, gridLoc='CC')
		q = sdiag(1/mesh.vol)*np.zeros(mesh.nC)
		q[txind] = 1.

		return q

class AcousticRx(Survey.BaseRx):

	def __init__(self, locs, **kwargs):
		self.locs = locs

	# Question: Why this does not work?
	# def getP(self):
	# 	print 'kang'
	# 	return mesh.getInterpolationMat(self.locs, 'CC')

	@property
	def nD(self):
		""" The number of data in the receiver."""
		return self.locs.shape[0]

	def getP(self, mesh):
		P = mesh.getInterpolationMat(self.locs, 'CC')
		return P



class SurveyAcoustic(Survey.BaseSurvey):
	"""
		**SurveyAcousitc**

		Geophysical Acoustic Wave data.

	"""

	def __init__(self, txList,**kwargs):
		self.txList = txList
		Survey.BaseSurvey.__init__(self, **kwargs)

	def projectFields(self, u):
		data = []

		for i, tx in enumerate(self.txList):
			Proj = tx.rxList[0].getP(self.prob.mesh)
			data.append(Proj*u[i])
		return data


class AcousticProblemSponge(Problem.BaseProblem):
	"""

	"""
	surveyPair = SurveyAcoustic
	Solver     = Solver
	storefield = True
	verbose = False
	stability = False
	sig = False

	def __init__(self, mesh, **kwargs):
		Problem.BaseProblem.__init__(self, mesh)
		self.mesh.setCellGradBC('dirichlet')
		Utils.setKwargs(self, **kwargs)

	def setSpongeBC(self, npad, dt, bcflag="all"):
		#TODO: performance of abosrbing
		self.bcflag = bcflag
		ax = self.mesh.vectorCCx[-npad]
		ay = self.mesh.vectorCCy[-npad]
		if bcflag == 'all':
			indy = np.logical_or(self.mesh.gridCC[:,1]<=-ay, self.mesh.gridCC[:,1]>=ay)
			indx = np.logical_or(self.mesh.gridCC[:,0]<=-ax, self.mesh.gridCC[:,0]>=ax)
		elif bcflag =='up':
			indy = self.mesh.gridCC[:,1]<=-ay
			indx = np.logical_or(self.mesh.gridCC[:,0]<=-ax, self.mesh.gridCC[:,0]>=ax)
		elif bcflag =='left':
			indy = self.mesh.gridCC[:,1]<=-ay
			indx = self.mesh.gridCC[:,0]>=ax
		elif bcflag =='right':
			indy = self.mesh.gridCC[:,1]<=-ay
			indx = self.mesh.gridCC[:,0]<=-ax
		else:
			raise Exception("Not implented!!")
		tempx = np.zeros_like(self.mesh.gridCC[:,0])
		tempx[indx] = (abs(self.mesh.gridCC[:,0][indx])-ax)**2
		tempx[indx] = tempx[indx]-tempx[indx].min()
		tempx[indx] = tempx[indx]/tempx[indx].max()
		tempy = np.zeros_like(self.mesh.gridCC[:,1])
		tempy[indy] = (abs(self.mesh.gridCC[:,1][indy])-ay)**2
		tempy[indy] = tempy[indy]-tempy[indy].min()
		tempy[indy] = tempy[indy]/tempy[indy].max()
		temp = tempx+tempy
		temp[temp>1.] = 1.

		f = 1.- temp*0.1
		self.sig = (1.-f)/f*2./dt


	def stabilitycheck(self, v, time, fmain):

		self.dxmin = min(self.mesh.hx.min(), self.mesh.hy.min())
		self.topt = self.dxmin/v.max()*0.5
		self.dt = time[1]-time[0]
		self.fmain = fmain
		self.wavelen = v.min()/self.fmain
		self.G = self.wavelen/self.dxmin

		if self.dt > self.topt:
			print "Warning: dt is greater than topt"
			self.stability = False
		elif self.G < 16.:
			print "Warning: Wavelength per cell (G) should be greater than 16"
			self.stability = False
		else:
			print "You are good to go:)"
			self.stability = True

		print ">> Stability information"
		print ("   dt: %5.2e s")%(self.dt)
		print ("   Optimal dt: %5.2e s")%(self.topt)
		print ("   Cell per wavelength (G): %5.2e")%(self.G)
		print ("   Optimal G: %5.2e")%(16.)

	def fields(self, v):

		if self.bcflag =='up':
			self.mesh.setCellGradBC(['dirichlet', ['dirichlet', 'neumann']])
		elif self.bcflag =='left':
			self.mesh.setCellGradBC([['neumann', 'dirichlet'], ['dirichlet', 'neumann']])
		elif self.bcflag =='right':
			self.mesh.setCellGradBC([['dirichlet', 'neumann'], ['dirichlet', 'neumann']])

		Grad = self.mesh.cellGrad
		Div = self.mesh.faceDiv
		AvF2CC = self.mesh.aveF2CC
		rho = 0.27*np.ones(self.mesh.nC)
		mu = rho*v**2
		DSig = sdiag(self.sig)
		Drhoi = sdiag(1/rho)
		MfmuiI = sdiag(1/(AvF2CC.T*(1/mu)))

		if self.stability==False:
			raise Exception("Stability condition is not satisfied!!")
		elif self.sig is False:
			print "Warning: Absorbing boundary condition was not set yet!!"
		start = clock()
		print ">> Start Computing Acoustic Wave"
		print (">> dt: %5.2e s")%(self.dt)
		print (">> Optimal dt: %5.2e s")%(self.topt)
		print (">> Main frequency, fmain: %5.2e Hz")%(self.fmain)
		print (">> Cell per wavelength (G): %5.2e")%(self.G)


		if self.storefield==True:
			P = []
			#TODO: parallize in terms of sources
			ntx = len(self.survey.txList)
			for itx, tx in enumerate(self.survey.txList):
				print ("  Tx at (%7.2f, %7.2f): %4i/%4i")%(tx.loc[0], tx.loc[0], itx+1, ntx)
				pn = np.zeros(self.mesh.nC)
				p0 = np.zeros_like(pn)
				un = np.zeros(self.mesh.nF)
				u0 = np.zeros_like(un)
				time = tx.time
				dt = tx.dt
				p = np.zeros((self.mesh.nC, time.size))
				q = tx.getq(self.mesh)
				for i in range(time.size-1):
					sn = tx.Wave(i+1)
					s0 = tx.Wave(i)
					pn = p0-dt*DSig*p0+Drhoi*dt*(Div*un+(sn-s0)/dt*q)
					p0 = pn.copy()
					# un = u0 - dt*sdiag(AvF2CC.T*self.sig)*u0 + dt*MfmuiI*Grad*p0
					un = u0 + dt*MfmuiI*Grad*p0
					u0 = un.copy()
					p[:,i+1] = pn

				P.append(p)
			elapsed = clock()-start
			print (">>Elapsed time: %5.2e s")%(elapsed)

			return P

		elif self.storefield==False:

			Data = []

			ntx = len(self.survey.txList)
			for itx, tx in enumerate(self.survey.txList):
				print ("  Tx at (%7.2f, %7.2f): %4i/%4i")%(tx.loc[0], tx.loc[0], itx+1, ntx)
				pn = np.zeros(self.mesh.nC)
				p0 = np.zeros_like(pn)
				un = np.zeros(self.mesh.nF)
				u0 = np.zeros_like(un)
				time = tx.time
				dt = tx.dt
				data = np.zeros((time.size, tx.nD))
				q = tx.getq(self.mesh)
				Proj = tx.rxList[0].getP(self.mesh)
				for i in range(time.size-1):
					sn = tx.Wave(i+1)
					s0 = tx.Wave(i)
					pn = p0-dt*DSig*p0+Drhoi*dt*(Div*un+(sn-s0)/dt*q)
					p0 = pn.copy()
					# un = u0 - dt*sdiag(AvF2CC.T*self.sig)*u0 + dt*MfmuiI*Grad*p0
					un = u0 + dt*MfmuiI*Grad*p0
					u0 = un.copy()
					data[i,:] =  Proj*pn

				Data.append(data)

			return Data

class AcousticProblemPML(Problem.BaseProblem):
	"""

	"""
	surveyPair = SurveyAcoustic
	Solver     = Solver
	storefield = True
	verbose = False
	stability = False
	sigx = False

	def __init__(self, mesh, **kwargs):
		Problem.BaseProblem.__init__(self, mesh)
		self.mesh.setCellGradBC('dirichlet')
		Utils.setKwargs(self, **kwargs)

	def setPMLBC(self, npad, dt, bcflag="all"):
		#TODO: performance of abosrbing
		self.bcflag = bcflag
		ax = self.mesh.vectorCCx[-npad]
		ay = self.mesh.vectorCCy[-npad]
		if bcflag == 'all':
			indy = np.logical_or(self.mesh.gridCC[:,1]<=-ay, self.mesh.gridCC[:,1]>=ay)
			indx = np.logical_or(self.mesh.gridCC[:,0]<=-ax, self.mesh.gridCC[:,0]>=ax)
		elif bcflag =='up':
			indy = self.mesh.gridCC[:,1]<=-ay
			indx = np.logical_or(self.mesh.gridCC[:,0]<=-ax, self.mesh.gridCC[:,0]>=ax)
		elif bcflag =='left':
			indy = self.mesh.gridCC[:,1]<=-ay
			indx = self.mesh.gridCC[:,0]>=ax
		elif bcflag =='right':
			indy = self.mesh.gridCC[:,1]<=-ay
			indx = self.mesh.gridCC[:,0]<=-ax
		else:
			raise Exception("Not implented!!")
		tempx = np.zeros_like(self.mesh.gridCC[:,0])
		tempx[indx] = (abs(self.mesh.gridCC[:,0][indx])-ax)**2
		tempx[indx] = tempx[indx]-tempx[indx].min()
		tempx[indx] = tempx[indx]/tempx[indx].max()
		tempy = np.zeros_like(self.mesh.gridCC[:,1])
		tempy[indy] = (abs(self.mesh.gridCC[:,1][indy])-ay)**2
		tempy[indy] = tempy[indy]-tempy[indy].min()
		tempy[indy] = tempy[indy]/tempy[indy].max()

		fx = 1.-tempx*0.1
		fy = 1.-tempy*0.1
		self.sigx = (1-fx)/fx*2./dt
		self.sigy = (1-fy)/fy*2./dt

	def stabilitycheck(self, v, time, fmain):

		self.dxmin = min(self.mesh.hx.min(), self.mesh.hy.min())
		self.topt = self.dxmin/v.max()*0.5
		self.dt = time[1]-time[0]
		self.fmain = fmain
		self.wavelen = v.min()/self.fmain
		self.G = self.wavelen/self.dxmin

		if self.dt > self.topt:
			print "Warning: dt is greater than topt"
			self.stability = False
		elif self.G < 16.:
			print "Warning: Wavelength per cell (G) should be greater than 16"
			self.stability = False
		else:
			print "You are good to go:)"
			self.stability = True

		print ">> Stability information"
		print ("   dt: %5.2e s")%(self.dt)
		print ("   Optimal dt: %5.2e s")%(self.topt)
		print ("   Cell per wavelength (G): %5.2e")%(self.G)
		print ("   Optimal G: %5.2e")%(16.)

	def fields(self, v):

		if self.bcflag =='up':
			self.mesh.setCellGradBC(['dirichlet', ['dirichlet', 'neumann']])
		elif self.bcflag =='left':
			self.mesh.setCellGradBC([['neumann', 'dirichlet'], ['dirichlet', 'neumann']])
		elif self.bcflag =='right':
			self.mesh.setCellGradBC([['dirichlet', 'neumann'], ['dirichlet', 'neumann']])

		Grad = self.mesh.cellGrad
		AvF2CC = self.mesh.aveF2CC
		AvF2CCv = self.mesh.aveF2CCV
		rho = 0.27*np.ones(self.mesh.nC)
		mu = rho*v**2
		Divvec = sp.block_diag((self.mesh.faceDivx, self.mesh.faceDivy))
		Mrhocc = sp.block_diag((sdiag(rho), sdiag(rho)))
		Mmuifvec = sdiag(AvF2CC.T*(1/mu))
		Msigf = sdiag(AvF2CCv.T*np.r_[self.sigx, self.sigy])
		Msigcc =  sp.block_diag((sdiag(self.sigx), sdiag(self.sigy)))
		MrhoccI = sdInv(Mrhocc)
		MmuifvecI = sdInv(Mmuifvec)
		Ivec = sp.hstack((sdiag(np.ones(self.mesh.nC)), sdiag(np.ones(self.mesh.nC))))

		if self.stability==False:
			raise Exception("Stability condition is not satisfied!!")
		elif self.sigx is False:
			print "Warning: Absorbing boundary condition was not set yet!!"
		start = clock()
		print ">> Start Computing Acoustic Wave"
		print (">> dt: %5.2e s")%(self.dt)
		print (">> Optimal dt: %5.2e s")%(self.topt)
		print (">> Main frequency, fmain: %5.2e Hz")%(self.fmain)
		print (">> Cell per wavelength (G): %5.2e")%(self.G)


		if self.storefield==True:
			Phi = []
			#TODO: parallize in terms of sources
			ntx = len(self.survey.txList)
			for itx, tx in enumerate(self.survey.txList):
				print ("  Tx at (%7.2f, %7.2f): %4i/%4i")%(tx.loc[0], tx.loc[0], itx+1, ntx)
				phi = np.zeros((self.mesh.nC, tx.time.size))
				phin = np.zeros(self.mesh.nC*2)
				phi0 = np.zeros_like(phin)
				un = np.zeros(self.mesh.nF)
				u0 = np.zeros_like(un)
				time = tx.time
				dt = tx.dt
				q = tx.getq(self.mesh)
				qvec = np.r_[q, q]*1/2

				for i in range(time.size-1):
					sn = tx.Wave(i+1)
					s0 = tx.Wave(i)
					phin = phi0-dt*(Msigcc*phi0)+dt*MrhoccI*(1/dt*(sn-s0)*qvec+Divvec*un)
					phi0 = phin.copy()
					un = u0 - dt*Msigf*u0 + dt*MmuifvecI*Grad*(Ivec*phi0)
					u0 = un.copy()
					phi[:,i+1] = Ivec*phi0
				Phi.append(phi)
			elapsed = clock()-start
			print (">>Elapsed time: %5.2e s")%(elapsed)

			return Phi

		elif self.storefield==False:

			Data = []

			ntx = len(self.survey.txList)
			for itx, tx in enumerate(self.survey.txList):
				print ("  Tx at (%7.2f, %7.2f): %4i/%4i")%(tx.loc[0], tx.loc[0], itx+1, ntx)
				phi = np.zeros((mesh.nC, time.size))
				phin = np.zeros(mesh.nC*2)
				phi0 = np.zeros_like(phin)
				un = np.zeros(mesh.nF)
				u0 = np.zeros_like(un)
				time = tx.time
				dt = tx.dt
				p = np.zeros((self.mesh.nC, time.size))
				q = tx.getq(self.mesh)
				qvec = np.r_[q, q]*1/2
				for i in range(time.size-1):
					sn = tx.Wave(i+1)
					s0 = tx.Wave(i)
					phin = phi0-dt*(Msigcc*phi0)+dt*MrhoccI*(1/dt*(sn-s0)*qvec+Divvec*un)
					phi0 = phin.copy()
					un = u0 - dt*Msigf*u0 + dt*MmuifvecI*Grad*(Ivec*phi0)
					u0 = un.copy()
					data[i,:] =  Proj*(Ivec*phi0)

				Data.append(data)

			elapsed = clock()-start
			print (">>Elapsed time: %5.2e s")%(elapsed)

			return Data

if __name__ == '__main__':

	time = np.linspace(0, 0.04, 2**9)
	dt = time[1]-time[0]
	options={'tlag':0.0025, 'fmain':400}
	rx = AcousticRx(np.vstack((np.r_[0, 1], np.r_[0, 1])))
	tx = AcousticTx(np.r_[0, 1], time, [rx], **options)
	survey = SurveyAcoustic([tx])
	wave = tx.RickerWavelet()
	cs = 0.5
	hx = np.ones(150)*cs
	hy = np.ones(150)*cs
	mesh = Mesh.TensorMesh([hx, hy], 'CC')
	prob = AcousticProblemPML(mesh)
	prob.pair(survey)
	prob.setPMLBC(30, dt, bcflag='all')
	prob.storefield = True
	v = np.ones(mesh.nC)*2000.
	prob.stabilitycheck(v, time, 100.)
	U = prob.fields(v)
	icount = 220
	extent = [mesh.vectorCCx.min(), mesh.vectorCCx.max(), mesh.vectorCCy.min(), mesh.vectorCCy.max()]
	fig, ax = plt.subplots(1, 1, figsize = (5,5))
	ax.imshow(np.flipud(U[0][:,icount].reshape((mesh.nCx, mesh.nCy), order = 'F').T), cmap = 'binary')
	plt.show()

