#!/usr/bin/python
# basis_functions.py
#
# Created by Travis Johnson on 2010-06-02.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *
from legendre_galerkin import *

N = 100;
x = linspace(0,2*pi, N)
figure(1)
for k in range(2):
	plot(x,real(exp(1j*k*x)),x,imag(exp(1j*k*x)))
	axis([0-.1, 2*pi+.1, -1.1, 1.1])
draw()
savefig('fourier.eps')
### Chebyshev
x=linspace(-1,1,N)
figure(2)
for k in range(4):
	plot(x,cos(k*arccos(x)))
	axis([min(x)-.1, max(x)+.1, -1.1, 1.1])
draw()
savefig('cheb.eps')
### Legendre
figure(3)
for k in range(4):
	y = zeros((N,1))
	for i in range(len(x)):
		y[i] = LegendrePolynomial(k,x[i])
	plot(x,y)
savefig('leg.eps')
### ModifiedLegendre
figure(4)
for k in range(4):
	y= zeros((N,1))
	for i in range(len(x)):
		y[i] = ModifiedLegendreBasis(k,x[i])
	plot(x,y)
savefig('mleg.eps')
N = 15
xUniform = linspace(-1,1,N+1)
xLG, w = LegendreGaussNodesAndWeights(N)
xLGL,w = LegendreGaussLobattoNodesAndWeights(N)
xCh, w = ChebyshevGaussNodesAndWeights(N)
xCGL,w = ChebyshevGaussLobattoNodesAndWeights(N)


figure(5),plot(xUniform,0*xUniform, 'o'),text( -.75, .025, "Uniform Points")
plot(xLG,  0*xLG+.1,'o'),text( -.75, .125, "Legendre-Gauss Points")
plot(xLGL,  0*xLGL+.2,'o'),text( -.75, .225, "Legendre-Gauss-Lobatto Points")
plot(xCh,  0*xCh+.3,'o'),text( -.75, .325, "Chebyshev-Gauss Points")
plot(xCGL,  0*xCGL+.4,'o'),text( -.75, .425,"Chebyshev-Gauss-Lobatto Points")
axis([-1, 1, -.025, .45])
savefig('quad_nodes.eps')
show()
