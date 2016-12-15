#!/usr/bin/python
# polycollocation.py
#
# Created by Travis Johnson on 2010-05-30.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *
import time
# algorithm 19: imported
# algorithm 25: imported

N, tMax = 12, .1;
smartTimeStep = int(floor( tMax/(2.51/(0.047 * N**4)) * 2))
print smartTimeStep
X,phi = LegendreCollocationIntegrator(N,smartTimeStep,4*N,.1, lambda x: sin(pi*(x+1)))
x = linspace(-1,1)
k = .80
exact = sin(pi*(x+1))*exp(-k**2*pi**2*.1)
plot(x,exact)
show()