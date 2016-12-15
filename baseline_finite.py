#!/usr/bin/python
# baseline_finite.py
#
# Created by Travis Johnson on 2010-06-01.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from scipy.linalg import *
import time



def AdvectionDiffusionTimeDerivative(phihat):
	L, nu = 2*pi, .2
	phihatdot = (nu*dot(D2x,phihat))
	return phihatdot

## algorithm 45
def FiniteStep(tn, dt, phihat):
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

def FiniteDriver(N, Nt, T,initialValues):
	dt = T/Nt
	print "timestep is %f"%(dt)
	tn = 0
	phi = initialValues
	for n in range(Nt+1):
		phinew = FiniteStep(tn, dt, phi)
		if n%(Nt//10) ==0:
			plot(phinew),title("time t=%f"%(tn)),draw()
		tn = (n+1)*dt
		phi = phinew
	return phi

nVec = 2**arange(2,10)
errvec = []
for N in nVec:
	x=linspace(0.,2.,N)
	dx = x[1]-x[0]
	firstrow, firstcol = zeros((N,1)), zeros((N,1))
	firstrow[1]= 1.
	firstcol[1]= -1.
	D1x =toeplitz(firstcol,firstrow)
	D1x[0,0],D1x[0,1],D1x[0,2] = -3,4, -1
	D1x[-1,-1],D1x[-1,-2],D1x[-1,-3] = 3,-4, 1
	D1x = D1x/(2*dx)
	firstrow, firstcol = zeros((N,1)), zeros((N,1))
	firstrow[0],firstrow[1] = -2.,1.
	firstcol[0],firstcol[1] = -2.,1.
	D2x = toeplitz(firstcol,firstrow)
	D2x[0,0], D2x[0,1], D2x[0,2], D2x[0,3] = 2./(1), -5./(1), 4./(1), -1./(1)
	D2x[-1,-1], D2x[-1,-2], D2x[-1,-3], D2x[-1,-4] = 2./(1),  -5./(1), 4./(1), -1./(1)
	D2x = D2x/dx**2
	t= .5
	evals = eigvals(D2x)
	emax = max(evals)
	print emax
	phi = FiniteDriver(N,20000,t, sin(pi*(x)))
	phi_exact = sin(pi*x)*exp(-.2*pi**2*t)
	errvec.append(linalg.norm(phi_exact - phi))
	print errvec

close()
figure(1),semilogy(nVec, array(errvec))
title('semilog error')
#legend(('err 1 deriv','err 2 deriv'))
figure(2),loglog(nVec, array(errvec))
title('loglog error of finite difference')
xlabel('N')
ylabel('$\log_{10} (\max \text{error})$')
savefig('finitediff_error.eps')
#legend(('err 1 deriv','err 2 deriv'))

show()


