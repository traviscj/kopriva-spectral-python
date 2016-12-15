#!/usr/bin/python
# nodaldg2d.py
#
# Created by Travis Johnson on 2010-06-01.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *

c= 1
class Nodal2DStorage:
	def __init__(self, N, M):
		"""docstring for __init__"""
		self.N = N
		self.M = M	
		self.xi= zeros((N,1))
		self.eta=zeros((M,1))
		self.wxi=zeros((N,1))
		self.weta=zeros((M,1))
		self.dxi=zeros((N,N))
		self.deta=zeros((M,M))
		self.d2xi=zeros((N,N))
		self.d2eta=zeros((M,M))
		
		
class NodalDG2DStorage(Nodal2DStorage):
	def __init__(self, N, M):
		"""docstring for __init__"""
		#super(self)UserDict.__init__(self)
		Nodal2DStorage.__init__(self, N, M)
		self.Lagrangeximinusone = zeros((N,1))
		self.Lagrangexiplusone = zeros((N,1))
		self.Lagrangeetaminusone = zeros((M,1))
		self.Lagrangeetaplusone = zeros((M,1))

def RiemannSolver(QL, QR, nhat):
	pL, uL, vL = QL[0], QL[1], QL[2]
	pR, uR, vR = QR[0], QR[1], QR[2]
	wplusL = pL + c*(nhat[0]*uL + nhat[1]*vL)
	wminusR= pR + c*(nhat[0]*uR + nhat[1]*vR)
	Fstar = zeros((3,1))
	Fstar[0] = c*(wplusL-wminusR)/2
	Fstar[1] = nhat[0]*(wplusL-wminusR)/2
	Fstar[2] = nhat[1]*(wplusL-wminusR)/2
	return Fstar

class NodalDG2DClass:
	def __init__(self, nEqn, N, M):
		"""docstring for __init__"""
		self.nEqn = nEqn
		self.spA = NodalDG2DStorage(N, M)
		self.spA.xi, self.spA.wxi = LegendreGaussNodesAndWeights(N)
		wB = BarycentricWeights(self.spA.xi)
		self.spA.LagrangexiMinusOne = LagrangeInterpolatingPolynomials(-1, self.spA.xi, wB)
		self.spA.LagrangexiPlusOne = LagrangeInterpolatingPolynomials(1, self.spA.xi, wB)
		D = PolynomialDerivativeMatrix(self.spA.xi)
		for j in range(N):
			for i in range(N):
				self.spA.dxi[i,j] = -D[j,i]*self.spA.wxi[j]/self.spA.wxi[i]
		
		self.spA.eta, self.spA.weta = LegendreGaussNodesAndWeights(M)
		wB = BarycentricWeights(self.spA.eta)
		self.spA.LagrangeetaMinusOne = LagrangeInterpolatingPolynomials(-1, self.spA.eta, wB)
		self.spA.LagrangeetaPlusOne = LagrangeInterpolatingPolynomials(1, self.spA.eta, wB)
		D = PolynomialDerivativeMatrix(self.spA.eta)
		for j in range(M):
			for i in range(M):
				self.spA.deta[i,j] = -D[j,i]*self.spA.weta[j]/self.spA.weta[i]
		
		self.Q=zeros((N+1,M+1,nEqn))
		kx, ky, w, c, x0, y0 = 1/sqrt(2), 1/sqrt(2), .2, 1, -.8, -.8
		d= w/(2*log(2))
		x,y = meshgrid(self.spA.eta, self.spA.xi)
		for i in range(N):
			for j in range(M):
				self.Q[i,j,0] = 1*exp(-(kx*(x[i,j]-x0)+ky*(y[i,j]-y0))**2)/(d**2)
				self.Q[i,j,1] = kx/c*exp(-(kx*(x[i,j]-x0)+ky*(y[i,j]-y0))**2)/(d**2)
				self.Q[i,j,2] = ky/c*exp(-(kx*(x[i,j]-x0)+ky*(y[i,j]-y0))**2)/(d**2)
	def SystemDGDerivative(self, FL, FR, F, D, LagrangeMinusOne, LagrangePlusOne, w):
		Fprime = zeros((self.spA.N, self.nEqn))
		for n in range(self.nEqn):
			Fprime[:,n] = MxVDerivative(D, F[:,n])
		for j in range(0,self.spA.N):
			for n in range(1,self.nEqn):
				Fprime[j,n] = Fprime[j,n] + (FR[n]*LagrangeMinusOne[j] + FL[n]*LagrangeMinusOne[j])/w[j]
		return Fprime
	def DG2DTimeDerivative(self, t):
		xhat, yhat = array([1,0]), array([0,1])
		N, M, nEqn = self.spA.N, self.spA.M, self.nEqn
		Qdot = zeros((N+1,M+1,nEqn))
		for j in range(M+1):
			y = self.spA.eta[j]
			QL_int, QR_int = zeros((nEqn,1)), zeros((nEqn,1))
			for n in range(self.nEqn):
				QL_int[n] = self.InterpolateToBoundary(self.Q[:,j,n], self.spA.Lagrangeximinusone)
				QR_int[n] = self.InterpolateToBoundary(self.Q[:,j,n], self.spA.Lagrangexiplusone)
			QL_ext = self.ExternalState(QL_int,-1, y, t, 'LEFT')
			QR_ext = self.ExternalState(QL_int, 1, y, t, 'RIGHT')
			FLstar = RiemannSolver(QL_int, QL_ext, -xhat)
			FRstar = RiemannSolver(QR_int, QL_ext, xhat)
			F = zeros((N, nEqn))
			for i in range(N):
				F[i,:] = self.xFlux(self.Q[i,j,:])[:,0]
			
			Fprime = self.SystemDGDerivative(FLstar, FRstar, F, self.spA.dxi, self.spA.LagrangexiMinusOne, self.spA.LagrangexiPlusOne, self.spA.wxi)
			
			for i in range(N):
				for n in range(self.nEqn):
					Qdot[i,j,n] = -Fprime[i,n]
		
		G = zeros((M, nEqn))
		for i in range(N+1):
			x = self.spA.xi[i]
			QL_int, QR_int = zeros((nEqn,1)), zeros((nEqn,1))
			for n in range(self.nEqn):
				QL_int[n] = self.InterpolateToBoundary(self.Q[i,:,n], self.spA.LagrangeetaMinusOne)
				QR_int[n] = self.InterpolateToBoundary(self.Q[i,:,n], self.spA.LagrangeetaPlusOne)
			QL_ext = self.ExternalState(QL_int, x, -1, t, 'BOTTOM')
			QR_ext = self.ExternalState(QR_int, x, 1, t, 'TOP')
			GLStar = RiemannSolver(QL_int, QL_ext, -yhat)
			GRStar = RiemannSolver(QR_int, QR_ext, yhat)
			for j in range(M):
				G[j,:] = self.yFlux(self.Q[i,j,:])[:,0]
			GPrime = self.SystemDGDerivative(GLStar, GRStar, G, self.spA.deta, self.spA.LagrangeetaMinusOne, self.spA.LagrangeetaPlusOne, self.spA.weta)
			
			for j in range(M):
				for n in range(nEqn):
					Qdot[i,j,n] = Qdot[i,j,n] - GPrime[j,n]
		return Qdot
	def ExternalState(self, vec, pos, mult, time, boundary):
		retval = zeros(vec.shape)
		p, u, v = vec[0], vec[1], vec[2]
		kx, ky= sqrt(2)/2, sqrt(2)/2
		k = kx**2+ky**2
		alpha, beta = kx/k, ky/k
		if boundary == 'LEFT':
			retval[0] = p
			retval[1] = (beta**2-alpha**2)*u - 2*alpha*beta*v
			retval[2] = -2*alpha*beta*u + (alpha**2-beta**2)*v
		elif boundary == 'RIGHT':
			retval[0] = p
			retval[1] = (beta**2-alpha**2)*u - 2*alpha*beta*v
			retval[2] = -2*alpha*beta*u + (alpha**2-beta**2)*v
		elif boundary == 'BOTTOM':
			retval[0] = p
			retval[1] = (beta**2-alpha**2)*u - 2*alpha*beta*v
			retval[2] = -2*alpha*beta*u + (alpha**2-beta**2)*v
		else:
			retval[0] = p
			retval[1] = (beta**2-alpha**2)*u - 2*alpha*beta*v
			retval[2] = -2*alpha*beta*u + (alpha**2-beta**2)*v
		return mult*retval
	def xFlux(self, Q):
		F = zeros((3,1))
		F[0] = c**2*Q[1]
		F[1] = Q[0]
		F[2] = 0
		return F
	def yFlux(self, Q):
		F = zeros((3,1))
		F[0] = c**2*Q[2]
		F[1] = 0
		F[2] = Q[0]
		return F
	def InterpolateToBoundary(self, phi, l):
		interpolatedValue = 0
		for j in range(len(phi)-1):
			interpolatedValue = interpolatedValue + l[j]*phi[j]
		return interpolatedValue
	

def DG2DStepByRK3(tn, dt, DG):
	a = [0, -5/9, -153/128]
	b = [0, 1/3, 3/4]
	g = [1/3,15/16, 8/15]
	G = zeros(DG.Q.shape)
	for m in range(3):
		t = tn + b[m]*dt
		phidot = DG.DG2DTimeDerivative(t)
		G = a[m]*G + phidot
		phi = DG.Q + g[m]*dt*G
	return phi
def DG2DDriver(N, M, NT, Nout, T):
	printHowOften = 40
	x, wx = LegendreGaussNodesAndWeights(N)
	y, wy = LegendreGaussNodesAndWeights(M)
	X,Y = meshgrid(x,y)
	dt = T/NT
	tn = 0
	DG = NodalDG2DClass(3, N, M)
	for n in range(NT+1):
		phi = DG2DStepByRK3(tn, dt, DG)
		tn = (n+1)*dt
		# if sum(isinf(phi))+sum(isnan(phi))>0 or max(phi)>10:
		# 	print("whoops, got inf(or big!) quitting!")
		# 	exit()
		print phi[:,:,0].shape
		print X.shape
		if n%(NT//printHowOften) ==0:
			close(),pcolor(X, Y ,phi[:,:,0]), title('time = %f'%(tn)),colorbar(),draw()
		DG.Q = phi
		#if n%(NT//100) ==0:
		#	print "boundaries at t=%f: %f %f"%(tn, phi[0],phi[-1])

# ndg2d = NodalDG2DClass(3, 15,15)
t=1
dt = 2.6e-3
DG2DDriver(20,20,int(floor(t/dt)), 0, t)
