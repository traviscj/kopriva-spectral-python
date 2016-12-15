#!/usr/bin/python
# nodal_discont_galerkin.py
#
# Created by Travis Johnson on 2010-06-01.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *

from KoprivaMethods import *



N, dt, t=36, 1.5e-4, 2

x, phi = DGDriver(N, int(floor(t/dt)), 0, t, initialValues)

plot(x,phi, x, initialValues(x-t))

	
# N, sigma = 20, .2
# x= linspace(-1,1,N)
# 
# f = cos(pi*x)
# xGL, wGL = LegendreGaussLobattoNodesAndWeights(N-1)
# w = BarycentricWeights(x)
# interpolant =zeros(len(x))
# for i in range(len(interpolant)):
# 	interpolant[i] = LagrangeInterpolation(xGL[i], x, f, w)
# print f
# print interpolant
# plot(x,f,xGL,interpolant,'x')
# show();