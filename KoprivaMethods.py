#!/usr/bin/python
# KoprivaMethods.py
#
# Created by Travis Johnson on 2010-05-30.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from scipy.special import legendre
import time

cot = lambda x: 1/tan(x)


## algorithm 18
def FourierDerivativeMatrix(N):
	D = zeros((N,N))
	for i in range(N):
		D[i,i] = 0
		for j in range(N):
			if j!=i:
				D[i,j]=.5*(-1)**(i+j)*cot((i-j)*pi/N)
				D[i,i]=D[i,i]-D[i,j]
	return D

## algorithm 19
def MxVDerivative(D,f):
	In = zeros(f.shape)
	for i in range(0,len(f)):
		t=0
		for j in range(0,len(f)):
			t += D[i,j]*f[j]
		In[i] = t
	return In

## algorithm 20
def LegendrePolynomial(k,x):
	if k==0: return 1
	if k==1: return x
	Lkm2 = 1
	Lkm1 = x
	for j in range(2, k+1):
		Lk = (2*j-1)/j*x*Lkm1 - (j-1)/j*Lkm2
		Lkm2 = Lkm1
		Lkm1 = Lk
	return Lk
## algo 21
def ChebyshevPolynomial(k,x):
	"""docstring for ChebyshevPolynomial"""
	Ks = 10
	if k==0:
		return 1
	if k==1:
		return x
	if k < Ks:
		Tkm2 = 1
		Tkm1 = x
		for j in range(2,k+1):
			Tk = 2*x*Tkm1 - Tkm2
			Tkm2 = Tkm1
			Tk = Tkm1
	else:
		Tk = cos(k*arccos(x))
	return Tk

## algo 22
def LegendrePolynomialAndDerivative(N, x):

	Ln, Lpn = 0,0
	if N==0:
		Ln, Lpn = 1,0
	elif N==1:
		Ln, Lpn = x,1
	else:
		Lnm2 = 1
		Lnm1 = x
		Lpnm2 = 0
		Lpnm1 = 1
		for k in range(2,N+1):
			Ln = (2*k-1)/k*x*Lnm1 - (k-1)/k*Lnm2
			Lpn = Lpnm2 + (2*k-1)*Lnm1
			Lnm2 = Lnm1
			Lnm1 = Ln
			Lpnm2 = Lpnm1
			Lpnm1 = Lpn
	return Ln, Lpn
## algo 23
def LegendreGaussNodesAndWeights(N):
	"""docstring for LegendreGaussNodesAndWeights"""
	nit, TOL = 50, 4*finfo(float).eps
	x = zeros((N+1,1))
	w = zeros((N+1,1))
	if N==0:
		x[0]=0
		w[0]=2
	elif N==1:
		x[0]= -1/sqrt(3)
		w[0] = 1
		x[1] = -x[0]
		w[1] = w[0]
	else:
		for j in range((N+1)//2):
			x[j] = -cos((2*j+1)/(2*N+2)*pi)
			for k in range(nit):
				Lnp1, Lpnp1 = LegendrePolynomialAndDerivative(N+1, x[j])
				
				delta = -Lnp1/Lpnp1
				x[j] = x[j] + delta
				if k > 5:
					print("something strange this way comes.")
				if abs(delta)<= TOL*abs(x[j]):
					break
			Lnp1, Lpnp1 = LegendrePolynomialAndDerivative(N+1, x[j])
			x[N-j] = -x[j]
			w[j] = 2/((1-x[j]**2)*(Lpnp1)**2)
			w[N-j] = w[j]
	
	if N%2 ==0:
		Lnp1, Lpnp1 = LegendrePolynomialAndDerivative(N+1, 0.0)
		x[N//2] = 0
		w[N//2] = 2/Lpnp1**2
	
	return x, w

## algorithm 25 - LegendreGaussLobattoNodesAndWeights
# calculate the points and weights for Gauss-Lobatto integration
def LegendreGaussLobattoNodesAndWeights(N):
	coeffs = array(legendre(N))
	deriv = arange(N,0,-1)
	load = []
	for elt in range(len(deriv)):
	    load.append(coeffs[elt]*deriv[elt])
	rts = roots(load)
	rts = concatenate((array([-1]), rts))
	rts = concatenate((rts,array([1])))
	x = sort(real(rts))
	w = 2/(N*(N+1)*polyval(coeffs, x)**2)
	return (x,w)
## algorithm 26
def ChebyshevGaussNodesAndWeights(N):
	x, w = zeros((N+1, 1)), zeros((N+1,1))
	for j in range(N+1):
		x[j] = -cos((2*j+1)/(2*N+2)*pi)
		w[j] = pi/(N+1)
	return x,w
## algorithm 27
def ChebyshevGaussLobattoNodesAndWeights(N):
	x, w = zeros((N+1, 1)), zeros((N+1,1))
	for j in range(N+1):
		x[j] = -cos(j/N*pi)
		w[j] = pi/(N)
	w[0] = w[0] / 2
	w[N] = w[N] / 2
	return x,w
## algorithm 30: 
def BarycentricWeights(x):
	N = len(x)
	w = ones((N,1))
	for j in range(1,N):
		for k in range(0,j):
			w[k] = w[k]*(x[k]-x[j])
			w[j] = w[j]*(x[j]-x[k])
			if w[k] == 0 or w[j] == 0:
				print "whoa, in barycentric!"
				exit()
	for j in range(N):
		w[j] = 1/w[j]
		if w[j] == 0:
			print "whoa, in barycentric!, %i %f"%(j, w[j])
			exit()
	return w
## algo 31
def LagrangeInterpolation(x,xj,f, w):
	numerator = 0
	denominator = 0
	for j in range(len(w)):
		if AlmostEqual(x,xj[j]):
			return f[j]
		t = w[j]/(x-xj[j])
		numerator += t*f[j]
		denominator += t
	return numerator/denominator

## algorithm 32:
def PolynomialInterpolationMatrix(x, w, xi):
	N = len(x)
	M = len(xi)
	T = zeros((M,N))
	for k in range(M):
		hasRowMatch = False
		for j in range(N):
			T[k,j] = 0
			if AlmostEqual(xi[k], x[j]):
				rowHasMatch = True
				T[k,j] = 1
		if rowHasMatch==False:
			s=0
			for j in range(N):
				t= w[j]/(xi[k]-x[j])
				T[k,j] = t
				s = s+t
			for j in range(N):
				T[k,j] = T[k,j]/s
	return T
## Algorithm 33:
def InterpolateToNewPoints(T,f):
	return dot(T,f)
## algorithm 37
def PolynomialDerivativeMatrix(x):
	N = len(x)
	D = zeros((N,N))
	w = BarycentricWeights(x)
	# print("barycentric weights: %s"%(str(w)))
	for i in range(0,N):
		D[i,i] = 0
		for j in range(N):
			if i != j:
				D[i,j] = w[j]/w[i]*1/(x[i]-x[j])
				D[i,i] -= D[i,j]
	return D

## algorithm 38
def mthOrderPolynomialDerivativeMatrix(m, x):
	N = len(x)
	w = BarycentricWeights(x)
	Dm = PolynomialDerivativeMatrix(x)
	if m==1:
		return Dm
	Dmm1 = Dm
	Dm = zeros((N,N))
	for k in range(2, m+1):
		print("taking %i to %i"%(k-1,k))
		for i in range(N):
			Dm[i,i] = 0
			for j in range(N):
				if i != j:
					Dm[i,j] = k/(x[i]-x[j])*(w[j]/w[i]*Dmm1[i,i] - Dmm1[i,j])
					Dm[i,i] = Dm[i,i] - Dm[i,j]
					if isnan(Dm[i,j]):
						print "A NaN! %i %i: "%(i,j)
						print "x: %f %f"%(x[i], x[j])
						print "w: %f %f"%(w[i], w[j])
						exit()
					elif isinf(Dm[i,j]):
						print "An Inf! %i %i: "%(i,j)
						print "x: %f %f"%(x[i], x[j])
						print "w: %f %f"%(w[i], w[j])
						exit()
		Dmm1 = Dm
	return Dm
## algorithm 41
def FourierCollocationTimeDerivative(phi, D):
	L, nu, N = 2*pi, .2, len(phi)
	k = fft.fftfreq(N, L/N/(2*pi))
	F = MxVDerivative(D, phi)
	F = nu*F - phi
	phidot = MxVDerivative(D,F)
	return phidot
## algorithm 42: renamed for naming consistency
def FourierCollocationStepByRK3(t_n,dt, phi, D):
	a=[0, -5/9, -153/128]
	b =[0, 1/3, 3/4]
	g =[1/3,15/16, 8/15]
	for m in range(3):
		t = t_n + b[m]*dt
		phidot = FourierCollocationTimeDerivative(phi,D)
		G = zeros(phi.shape)
		for j in range(len(phi)):
			G[j] = a[m]*G[j] + phidot[j]
			phi[j]=phi[j] + g[m]*dt*G[j]
	return phi

##  algorithm 43
# extra argument: initialValues(N)
def FourierCollocationDriver(N, Nt, T,initialValues):
	D = FourierDerivativeMatrix(N)
	dt = T/Nt
	tn = 0
	phi = initialValues(N)
	x=linspace(0,2*pi,N)
	for n in range(Nt+1):
		phinew = FourierCollocationStepByRK3(tn, dt, phi, D)
		#plot(phinew),title("time t=%f"%(tn)),draw()#, time.sleep(3)
		tn = (n+1)*dt
		phi = phinew
	return x, phi
	
## Algorithm 44
def AdvectionDiffusionTimeDerivative(phihat):
	L, nu = 2*pi, .2
	N = len(phihat)
	k = fft.fftfreq(N, L/N/(2*pi))
	phihatdot = -1*(1j*k+nu*k**2)*phihat
	return phihatdot

## algorithm 45
def FourierGalerkinStep(tn, dt, phihat):
	a = [0, -5/9, -153/128]
	b = [0, 1/3, 3/4]
	g = [1/3,15/16, 8/15]
	G = zeros(phihat.shape)
	for m in range(3):
		t = tn + b[m]*dt
		phihatdot = AdvectionDiffusionTimeDerivative(phihat)
		G = a[m]*G + phihatdot
		phihat = phihat + g[m]*dt*G
	return phihat

## algorithm 46: used fft.ifft
def EvaluateFourierGalerkinSolution(x,phihat):
	N = len(phihat)
	M = len(x)
	phi = zeros((M,1))
	for k in range(-N//2, N//2):
		phi[:,0] = phi[:,0] + real(phihat[k]/N*exp(1j*k*x))
	return phi

## algorithm 47
def FourierGalerkinDriver(N, NT, T, Nout,initialValues):
	dt = T/NT
	print("run with T=%f, NT=%i, dt=%f"%(T,NT,dt))
	tn = 0
	phihat = fft.fft(initialValues(N))
	x2=linspace(0,2*pi,N+1)
	x = x2[:-1]
	k = fft.fftfreq(N, 2*pi/N/(2*pi))
	for n in range(0,NT):
		phihat = FourierGalerkinStep(tn, dt, phihat)
		tn = (n+1)*dt
		if n%(NT//10) ==0:
			figure(1),plot(k,phihat,'.'),title('k-space t=%f'%tn),draw()
			figure(2),plot(x/pi,fft.ifft(phihat)),title('t = %f'%tn),draw()
			time.sleep(.75)
	dx = 2*pi/Nout
	x = linspace(0, 2*pi, Nout+1)
	#phi = fft.ifft(phihat)
	phi = EvaluateFourierGalerkinSolution(x,phihat)
	return (x,phi)

## algorithm 50
def CollocationStepByRK3(tn, dt, phi, D, timeDeriv, gL, gR):
	a=[0, -5/9, -153/128]
	b =[0, 1/3, 3/4]
	g =[1/3,15/16, 8/15]
	for m in range(3):
		t = tn + b[m]*dt
		phidot = timeDeriv(phi,D)
		G = zeros(phi.shape)
		for j in range(len(phi)):
			G[j] = a[m]*G[j] + phidot[j]
			phi[j]=phi[j] + g[m]*dt*G[j]
	phi[0] = gL
	phi[-1]= gR
	return phi
## algorithm 51:
def LegendreCollocationIntegrator(N, NT, Nout, T, initialValues):
	def TDerivative(phi, D):
		k = .8
		return k**2*dot(D,phi)
	x, w = LegendreGaussLobattoNodesAndWeights(N)
	D2 = mthOrderPolynomialDerivativeMatrix(2, x)
	dt = T/NT
	tn = 0
	phi = initialValues(x)
	vals = linalg.eigvals(D2)
	print max(vals)
	for n in range(NT+1):
		phi = CollocationStepByRK3(tn, dt, phi, D2, TDerivative, 0, 0)
		tn = (n+1)*dt
		if sum(isinf(phi))+sum(isnan(phi))>0 or max(phi)>100:
			print("whoops, got inf(or big!) quitting!")
			exit()
		plot(x, phi,'o'), title('time = %f'%(tn)),draw()
		time.sleep(.5);
	
	X=zeros((Nout,1))
	for j in range(Nout):
		X[j] = -1 + 2*j/Nout
	wBary = BarycentricWeights(x)
	T = PolynomialInterpolationMatrix(x,wBary,X)
	phi_interp = InterpolateToNewPoints(T,phi)
	return X,phi_interp

def LagrangeInterpolatingPolynomials(x, xj, w):
	xMatchesNode = False
	l = zeros((len(xj),1))
	for j in range(len(xj)):
		l[j] = 0.0
		if AlmostEqual(x,xj[j]):
			l[j] = 1.0
			xMatchesNode = True
	if xMatchesNode:
		return l
	s=0
	for j in range(len(xj)):
		t = w[j]/(x-xj[j])
		l[j] = t
		s = s+t
	for j in range(len(xj)):
		l[j] = l[j]/s
	return l

class NodalDiscontinuousGalerkin:
	def __init__(self, N, c, initialValues):
		self.N, self.c = N, c
		self.x, self.weights = LegendreGaussNodesAndWeights(N)
		#self.phi = zeros((N,N))
		self.phi = initialValues(self.x)
		wbary = BarycentricWeights(self.x)
		self.LagrangeMinusOne = LagrangeInterpolatingPolynomials(-1, self.x, wbary)
		self.LagrangePlusOne = LagrangeInterpolatingPolynomials(1, self.x, wbary)
		D = PolynomialDerivativeMatrix(self.x)
		self.Dhat = zeros((N+1,N+1))
		for j in range(self.N+1):
			for i in range(self.N+1):
				self.Dhat[i,j] = -D[j,i]*self.weights[j]/self.weights[i]
	def DGDerivative(self,phiL, phiR, phi):
		phiprime = MxVDerivative(self.Dhat, phi)
		for j in range(self.N+1):
			phiprime[j] = phiprime[j] + (phiR*self.LagrangePlusOne[j] - phiL*self.LagrangeMinusOne[j])/self.weights[j]
		return phiprime
	def InterpolateToBoundary(self, phi, l):
		interpolatedValue = 0
		for j in range(self.N+1):
			interpolatedValue = interpolatedValue + l[j]*phi[j]
		return interpolatedValue
	def DGTimeDerivative(self, t):
		def g(t):
			sigma = .2
			xi = 0
			#return 1
			return exp(-log(2)*(xi-t)**2/sigma**2)
		if self.c >0:
			phiL = g(t)
			phiR = self.InterpolateToBoundary(self.phi, self.LagrangePlusOne)
		else:
			phiR = g(t)
			phiL = self.InterpolateToBoundary(self.phi, self.LagrangeMinusOne)
		phidot = -self.c*self.DGDerivative(phiL, phiR, self.phi)
		return phidot

def DGStepByRK3(tn, dt, DG):
	a = [0, -5/9, -153/128]
	b = [0, 1/3, 3/4]
	g = [1/3,15/16, 8/15]
	G = zeros(DG.phi.shape)
	for m in range(3):
		t = tn + b[m]*dt
		phidot = DG.DGTimeDerivative(t)
		G = a[m]*G + phidot
		phi = DG.phi + g[m]*dt*G
	return phi
def DGDriver(N, NT, Nout, T, initialValues):
	x, w = LegendreGaussNodesAndWeights(N)
	dt = T/NT
	tn = 0
	DG = NodalDiscontinuousGalerkin(N, pi,initialValues)
	for n in range(NT+1):
		phi = DGStepByRK3(tn, dt, DG)
		tn = (n+1)*dt
		if sum(isinf(phi))+sum(isnan(phi))>0 or max(phi)>10:
			print("whoops, got inf(or big!) quitting!")
			exit()
		if n%(NT//10) ==0:
			close(),plot(x, phi,'.-'), title('time = %f'%(tn)),draw()
			savefig('DG_plot_tn%i.eps'%(n+1))
		DG.phi = phi
		if n%(NT//100) ==0:
			print "boundaries at t=%f: %f %f"%(tn, phi[0],phi[-1])

	# X=zeros((Nout,1))
	# 	for j in range(Nout):
	# 		X[j] = -1 + 2*j/Nout
	# 	wBary = BarycentricWeights(x)
	# 	T = PolynomialInterpolationMatrix(x,wBary,X)
	# 	phi_interp = InterpolateToNewPoints(T,phi)
	return x,phi
def initialValues(x):
	sigma = 0.2
	return exp(-log(2)*(x+1)**2/sigma**2)
## algorithm 139
def AlmostEqual(a,b):
	epsilon = finfo(float).eps
	if a==0 or b==0:
		if abs(a-b) <= 2*epsilon:
			return True
		else:
			return False
	else:
		if abs(a-b) <= epsilon*abs(a) and abs(a-b)<= epsilon*abs(b):
			return True
		else:
			return False
			
## algorithm 140
def TriDiagonalSolve(L, D, U, y):
	dhat = zeros((N,1))
	for j in range(N):
		dhat[j] = d[j]
	for j in range(1,N):
		dhat[j] = dhat[j] - L[j]/dhat[j-1]*U[j-1]
		y[j] = y[j] - L[j]/dhat[j-1]*y[j-1]
	x = zeros((N,1))
	x[N] = y[N]/hatd[N]
	for j in range(N-1,-1,-1):
		x[j] = (y[j]-u[j]*x[j+1])/d[j]
	return x