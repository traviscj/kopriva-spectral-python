#!/usr/bin/python
# fourier_collocation.py
#
# Created by Travis Johnson on 2010-05-30.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
import time

from KoprivaMethods import *

L = 2*pi

def initialValues(N):
	x = linspace(0,L-L/N,N)
	#return 3/(5-4*cos(x))
	return sin(x-pi)
	#return exp(-3*(x-L/2)**2)
	
nVec=2**arange(3,7)
for n in nVec:
	t, dt = 1, 2.5e-3
	
	x, phi = FourierCollocationDriver(n, int(floor(1/dt)), t, initialValues)
	phi_exact = sin(pi*(x+1))*exp(-.2*pi**2*1)
	k = fft.fftfreq(n, L/n/(2*pi))
	#phi_exact = fft.ifft(exp(-(1j*k+.2*k**2)*t*2)*fft.fft(initialValues(n)))
	close()
	#phi_exact = exp(-3*(x-L/2)**2)*exp(-.2*t)
	plot(x,phi,x,phi_exact,'o-')
	legend(('calc','exact'))
	draw()
	err = linalg.norm(phi-phi_exact)
	time.sleep(3)
	print err
