#!/usr/bin/python
# testmore.py
#
# Created by Travis Johnson on 2010-05-28.
# Copyright (c) 2010 . All rights reserved.

from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *

x, w  = LegendreGaussNodesAndWeights(6)
for i in range(len(x)):
	print "%9.16e\t %9.16e"%(x[i],w[i])