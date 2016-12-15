#!/usr/bin/python
# eigenvalues.py
#
# Created by Travis Johnson on 2010-06-02.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *

errlist = []
nVec=arange(32,32)
nVec=[31]
for N in nVec:
	xlgl, w = LegendreGaussNodesAndWeights(N)
	xcheb, w= ChebyshevGaussNodesAndWeights(N)
	D1cheb = mthOrderPolynomialDerivativeMatrix(2,xcheb)
	D1 = mthOrderPolynomialDerivativeMatrix(2,xlgl)
	eigs = eigvals(D1)
	eigscheb = eigvals(D1cheb)
	figure(1),plot(real(eigs),imag(eigs),'o')
	figure(2),plot(real(eigscheb),imag(eigscheb),'o')
show()