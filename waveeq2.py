#!/usr/bin/env python
# encoding: utf-8
"""
waveeq2.py

Created by Travis Johnson on 2010-05-23.
Copyright (c) 2010 . All rights reserved.
"""

import sys
import os
from chebdiff import cheb
from numpy import *
from pylab import plot, figure, legend, show

def main():
	N = 30
	L = 1
	print("init D1")
	(D1,x1)=cheb(N,-L,0)
	print("done")
	print("init D2")
	(D2,x2)=cheb(N,0,L)
	print("done")
	D = zeros((2*N+1,2*N+1))
	print D2.shape
	print D1.shape
	print D.shape
	D[N:,N:] = D2[:,:]
	D[:N+1,:N+1] = D1[:,:]
	x=hstack((hstack((x1[:-1],array([x1[-2]/10, 0, x2[1]/10]))),x2[1:]))
	x=hstack((hstack((x1[:-1],array([x1[-2], 0, x2[1]]))),x2[1:]))
	x=hstack((hstack((x1[:-1],array([ 0]))),x2[1:]))
	print x
	u=exp(x)*sin(5*x)
	print("will do mult...")
	ux=dot(D,u)
	print("done")
	ux_exact = exp(x)*(sin(5*x)+5*cos(5*x))
	print("plotting")
	plot(x,ux,x,ux_exact)
	legend(('calc','exact'))
	show()
	
	

if __name__ == '__main__':
	main()

