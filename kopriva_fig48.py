#!/usr/bin/python
# testmore.py
#
# Created by Travis Johnson on 2010-05-28.
# Copyright (c) 2010 . All rights reserved.

from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *

errlist = []
nVec=range(4,36,2)
for N in nVec:
	x, w = LegendreGaussLobattoNodesAndWeights(N)
	D2 = mthOrderPolynomialDerivativeMatrix(2, x)
	
	plot(x,-pi**2*sin(pi*x),'.',x,dot(D2,sin(pi*x)),'o')
	legend(('exact','approx'))
	err = max(abs(-pi**2*sin(pi*x)-dot(D2,sin(pi*x))))
	errlist.append(err)
	title('err = %e'%err)
	draw()
	time.sleep(1)

figure(2),semilogy(nVec,array(errlist))
show()