#!/usr/bin/python
# continuousgalerkin.py
#
# Created by Travis Johnson on 2010-05-31.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *

## algo 57
def CGDerivativeMatrix(N):
	x, w = LegendreGaussLobattoNodesAndWeights(N)
	D = PolynomialDerivativeMatrix(x)
	G = zeros((len(x),len(x)))
	for j in range(len(x)):
		for n in range(len(x)):
			s = 0
			for k in range(len(x)):
				s += D[k,n]*D[k,j] * w[k]
			G[j,n] = s/w[j]
	return G,x
## algorithm 50
def GalerkinStepByRK3(tn, dt, phi, D, timeDeriv):
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
	return phi
## algo 57a
def CGDerivativeMatrixIntegrator(N, NT, Nout, T, initialValues):
	def TDerivative(phi, D):
		k = .95
		df = -k**2*dot(D,phi)
		return df
	D,x =CGDerivativeMatrix(N)
	D[0,:], D[-1,:] = 0, 0
	#D[:,0], D[:,-1] = 0, 0
	print x
	dt = T/NT
	tn = 0
	phi = initialValues(x)
	vals = linalg.eigvals(D)
	print max(vals)
	for n in range(NT+1):
		phi = GalerkinStepByRK3(tn, dt, phi, D, TDerivative )
		print phi[0], phi[-1]
		tn = (n+1)*dt
		if sum(isinf(phi))+sum(isnan(phi))>0 or max(phi)>100:
			print("whoops, got inf(or big!) quitting!")
			exit()
		if n%(NT//10) ==0:
			exact = sin(pi*(x+1))*exp(-1/.95**2*pi**2*tn)
			diff = linalg.norm(exact - phi)
			plot(x, phi,'o',x,exact), title('time = %f, norm = %f'%(tn, diff)),draw()
			time.sleep(.5)
	
		
	X=zeros((Nout,1))
	for j in range(Nout):
		X[j] = -1 + 2*j/Nout
	wBary = BarycentricWeights(x)
	T = PolynomialInterpolationMatrix(x,wBary,X)
	phi_interp = InterpolateToNewPoints(T,phi)
	return X,phi_interp

N=12
G,x = CGDerivativeMatrix(N)
G = G[1:,1:]
lambdamax = max(eigvals(G))
print lambdamax
t, dt = .1, 1e-3
print "dt =",dt
NT = 380
x, phi = CGDerivativeMatrixIntegrator(25, NT, 25, t, lambda x: sin(pi*(x+1)))

# print len(x)
# print len(phi)
plot(x,phi)
savefig('continuousgalerkin.png')
