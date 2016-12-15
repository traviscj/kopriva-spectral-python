#!/usr/bin/env python
# encoding: utf-8
"""
deriv_test.py

Created by Travis Johnson on 2010-05-23.
Copyright (c) 2010 . All rights reserved.
"""

import sys
import os
from chebdiff import cheb
from numpy import *
from pylab import plot, figure, legend, show,savefig,title,semilogy

def main():
	# initialize main problem variables
	(N, L) = (20,1)
	# initialize chebyshev-spaced points, appropriate to the $i^{th}$ element
	(D1,x1)=cheb(N,-L,0)
	(D2,x2)=cheb(N,0,L)
	# build the big domain for error checking/etc
	x=hstack((hstack((x1[:],array([0]))),x2[:]))
	# starting with a nonperiodic function
	u=exp(x)*sin(5*x)
	# whose derivative we can calculate
	ux_exact = exp(x)*(sin(5*x)+5*cos(5*x))
	
	# build the main differentation matrix
	D=zeros((2*(N+1)+1,2*(N+1)+1));
	# by putting D1 in the upper left
	D[:N+1,:N+1] = D1[:,:]
	# and D2 in the lower right
	D[N+2:,N+2:] = D2[:,:]
	# this doesn't really matter, but it doesn't make sense to have all rows &
	# cols of this empty, and I think this makes sense for functions where
	# u(x=0) != 0
	D[N+1,N+1]=1;
	
	# calculate Du, which is u_x
	ux=dot(D,u)
	# correct the middle point by taking the average.
	ux[N+1] = .5*ux[N]+.5*ux[N+2]
	
	# plot the solution
	figure(1);plot(x,ux,x,ux_exact)
	legend(('calc','exact'))
	title('spectral element deriv approximation to u')
	savefig('approx_specelt_deriv.pdf')
	# and plot the error
	figure(2);semilogy(x,abs(ux-ux_exact))
	title('error of spectral element deriv approximation to u')
	savefig('approx_specelt_deriv_err.pdf')
	show()
	
	

if __name__ == '__main__':
	main()

