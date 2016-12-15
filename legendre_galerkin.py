#!/usr/bin/python
# legendre_galerkin.py
#
# Created by Travis Johnson on 2010-05-31.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *



## algo 52: ModifiedLegendreBasis
def ModifiedLegendreBasis(k,x):
	phik = (LegendrePolynomial(k,x) - LegendrePolynomial(k+2,x))/sqrt(4*k+6)
	return phik
	
## algo 53: EvaluateLegendreGalerkinSolution
def EvaluateLegendreGalerkinSolution(N,x,phihat):
	phi = 0
	for k in range(N-2):
		phi += phihat*ModifiedLegendreBasis(k,x)
	return phi

alpha = lambda n: 1/sqrt(4*n+6)
gamma = lambda n: -2/(2*n+1)
mu = lambda n: -2/(2*n+5)
beta = lambda n: -(gamma(n)+mu(n))
## algo 54: initTMatrix
def initTMatrix(n, p):
	diag = zeros((N,1))
	lower = zeros((N,1))
	upper = zeros((N,1))
	for j in range(N):
		diag[j] = beta(2*j+p)*alpha(2*j+p)**2
	for j in range(1,N):
		lower[j] = gamma(2*j+p)*alpha(2*j+p)*alpha(2*(j-1)+p)
		upper[j-1]=l[j]
	return (lower, diag, upper)

## algo 55: 
def ModifiedCoefsFromLegendreCoeffs(phihat):
	N = len(phihat)
	
	# even coeffs 
	M = int(floor((N-2)/2))
	(L, D, U) = initTMatrix(M, 0)
	rhs = zeros((M,1))
	for j in range(M):
		rhs[j] = mu(2*j)*alpha(2*j)*phihat(2*j+2) - alpha(2*j) * gamma(2*j) * phihat(2*j)
	b = TriDiagonalSolve(M, L, D, U, rhs)
	for j in range(M):
		phihat[2*j] = b[j]
	
	# odd coeffs
	M = int(floor((N-2+1)/2)-1)
	(L, D, U) = initTMatrix(M,1)
	for j in range(M):
		rhs[j] = mu(2*j+1)*alpha(2*j+1)*phihat(2*j+3) - alpha(2*j+1)*gamma(2*j+1)*phihat(2*j+1)
	b = TriDiagonalSolve(M, L, D, U, rhs)
	for j in range(M):
		phihat[2*j+1] = b[j]
	
	return phihat

# algo 56: LegendreGalerkinStep
def LegendreGalerkinStep(tn, dt, phihatn):
	N = len(phihat)
	
	# even coeffs
	M = int(floor(N-2)/2)
	L, D, U = initTMatrix(M,0)
	rhs = zeros((M,1))
	rhs[0] = (D[0] - dt/2)*phihatn[0] + U[0]*phihatn[2]
	for j in range(1,M-1):
		rhs[j] = L[j]*phihatn[2*(j-1)] + (D[j]-dt/2)*phihatn[2*j] + U[j]*phihatn[2*(j+1)]
	rhs[M] = (D[m]-dt/2)*phihatn[2*M] + L[M]*phihatn[2*(M-1)]
	for j in range(M):
		D[j] = D[j] + dt/2
	phihat = TriDiagonalSolve(M, L, D,U, rhs)
	for j in range(M):
		phihatnew[2*j] = phihat[j]
	
	# even coeffs
	M = int(floor(N-2+1)/2)
	L, D, U = initTMatrix(M,1)
	rhs = zeros((M,1))
	rhs[0] = (D[0] - dt/2)*phihatn[1] + U[0]*phihatn[3]
	for j in range(1,M-1):
		rhs[j] = L[j]*phihatn[2*(j-1)+1] + (D[j]-dt/2)*phihatn[2*j+1] + U[j]*phihatn[2*(j+1)+1]
	rhs[M] = (D[m]-dt/2)*phihatn[2*M+1] + L[M]*phihatn[2*(M-1)+1]
	for j in range(M):
		D[j] = D[j] + dt/2
	phihat = TriDiagonalSolve(M, L, D,U, rhs)
	for j in range(M):
		phihatnew[2*j+1] = phihat[j]

	return phihatnew
	

## algo unnumbered(post 56)
def LegendreGalerkinDriver(N,NT, T, Nout,initialValues):
	dt = T/NT
	tn = 0
	dx = 2*pi/Nout
	x = linspace(0, 2*pi, Nout+1)
	
	phihat = ModifiedCoefsFromLegendreCoefs(LegendreCoefs(initialvalues))
	for n in range(0,NT):
		phihat = LegendreGalerkinStep(tn, dt, phihat)
		if n%600 ==0:
			plot(EvaluateLegendreGalerkinSolution(N, x, phihat)),draw()
		tn = (n+1)*dt
	phi = EvaluateLegendreGalerkinSolution(N, x, phihat)
	return (x,phi)	
