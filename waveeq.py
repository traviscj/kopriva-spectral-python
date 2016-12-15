#!/usr/bin/env python
# encoding: utf-8
"""
waveeq.py

Created by Travis Johnson on 2010-05-23.
Copyright (c) 2010 . All rights reserved.
"""

import sys, os, time
from chebdiff import cheb
from numpy import *
from pylab import plot, figure, legend, show,savefig,title,semilogy,close,draw,axis

def main():
	errs=[]
	runForN(13,1)
def runForN(N, printTimes=0):
	# initialize main problem variables
	(L, dt, tMax, nImgs) = (1, 0.001, 5.0, 50)
	# initialize chebyshev-spaced points, appropriate to the $i^{th}$ element
	(D1,x1)=cheb(N,-L,0)
	(D2,x2)=cheb(N,0,L)
	# build the big domain for error checking/etc
	x=hstack((hstack((x1[:],array([0]))),x2[:]))
	# initial condition
	u0=zeros((len(x),1))
	u1=zeros((len(x),1))
	u1[0]=sin(pi*dt)
	
	# build the main differentation matrix
	D=zeros((2*(N+1)+1,2*(N+1)+1));
	# by putting D1 in the upper left
	D[:N+1,:N+1] = D1[:,:]
	# and D2 in the lower right
	D[N+2:,N+2:] = D2[:,:]
	# this doesn't really matter, but it doesn't make sense to have all rows &
	# cols of this empty, and I think this makes sense for functions where
	# u(x=0) != 0
	#D[N+1,N+1]=1;
	
	for n in range(2,int(floor(tMax/dt))):
		# precompute derivatives
		(Duxn,Duxnm1) = (dot(D,u1),dot(D,u0))
		# average at midpoint.
		(Duxn[N+1],Duxnm1[N+1]) = (1.*Duxn[N]+0*Duxn[N+2], 1.*Duxnm1[N]+0*Duxnm1[N+2])
		(Duxn[N+2],Duxnm1[N+2]) = (1.*Duxn[N+1],1.*Duxnm1[N+1])

		# evaluate newest time step
		u2 = u1 + dt*(-1.5*Duxn + .5*Duxnm1)
		# fix boundary condition
		u2[0] = sin(pi*dt*n)
		#u2[N+0] = u2[N+1]
		
		# slide n+1->n
		(u0,u1)=(u1,u2)
		
		# print a pretty picture, if appropriate
		if printTimes>0 and mod(n,int(tMax/dt/nImgs))==0:
			close(1);
			figure(1),plot(x,u2,'.-',x,sin(pi*(n*dt-x-1))*(x<= n*dt-1))
            # savefig('waveeq_sol_%0.4i.pdf'%(n));
            # close(1);
			axis([-1,1,-1,1])
			savefig('waveeq_sol_%0.4i.png'%(n));
			draw()
			time.sleep(1)
	return linalg.norm(u2.T-sin(pi*(n*dt-x-1))*(x<= n*dt-1),inf)

if __name__ == '__main__':
	main()

