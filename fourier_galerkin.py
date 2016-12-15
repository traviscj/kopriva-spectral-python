#!/usr/bin/python
# fourier_galerkin.py
#
# Created by Travis Johnson on 2010-05-30.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *

L=2*pi
def initialValues(N):
	x = linspace(0,L,N)
	return exp(-3*(x-L/2)**2)
	#return 3/(5-4*cos(x))

t,nu=2, .2
(x, phi) = FourierGalerkinDriver(15, 200, t, 50, initialValues)
print x
figure(3),plot(x/(pi),(phi),x/2,exp(-3*(x-L/2)**2)*exp(-nu*t))
legend(('calc','exact'))
show()
# y0 = initialValues(16)
# k = linspace(-len(y0)//2, len(y0)//2, len(y0))
# haty0 = fft.fft(y0)
# print k.shape
# print haty0.shape
# figure(2),plot(k,haty0,k,2**(-abs(k)))
# show()