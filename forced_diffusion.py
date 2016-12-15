#!/usr/bin/python
# forced_diffusion.py
#
# Created by Travis Johnson on 2010-05-28.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *
# L - Domain
# K - Number of Elements
# N - nodes per element
L, K, N = 16, 3, 10
TOL = 10**(-10)

#a=linspace(0, L, K+1)
a= array([-8, -3, 3, 8])
dx=zeros((K+1,1))
for i in range(1,len(dx)):
	print i
	dx[i] = (a[i]-a[i-1])
print a
print dx

def Gmat(j,m,w):
	summand = 0
	for l in range(N+1):
		Lj, Lpj = LegendrePolynomialAndDerivative(j,xi[l])
		Lm, Lpm = LegendrePolynomialAndDerivative(m,xi[l])
		summand += w[l]*Lpj*Lpm	
	return summand
	
xi, w= LegendreGaussLobattoNodesAndWeights(N)
G=zeros((N,N))
for i in range(N):
	for j in range(N):
		G[i,j] = Gmat(i,j,w)
globalG = zeros((K*N+1,K*N+1))
for k in range(K,0,-1):
	for i in range(N+1):
		for j in range(N+1):
			globalG[N*(k-1)+i,N*(k-1)+j] = Gmat(i,j,w)

globalG[-1,:], globalG[:,-1]=0,0
globalG[-1,-1]=1
globalG[0,0]=1



figure(1),spy(globalG)
x = zeros((K*(N)+1,1))
for k in range(1,K+1):
	for i in range(N+1):
		print k, i
		x[i+(k-1)*N] = a[k-1] + (xi[i]+1)/2*dx[k]
print x

#figure(2),plot(x,0*x,'o',array([a,a]),array([-1,1]))
#axis([-.1+min(a),.1+max(a),-1,1])

f = cos(pi/4*x)
print "the shape of f is ", f.shape
print "the shape of globalG is", globalG.shape
phi = dot(globalG,f)
df = zeros(phi.shape)
for i in range(len(df)):
	df = df +  phi*(LegendrePolynomialAndDerivative(i,xi)[1])

figure(3),plot(x,f,x,df,'o')
print min(x)
print max(x)
print min(f)
print max(f)
#axis((min(x), max(x), min(f), max(f)))
show()