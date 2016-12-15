#!/usr/bin/python
from __future__ import division
from numpy import *
from time import sleep
from pylab import plot,show,figure,semilogy,legend

# solve 
#	PDE		-u_{xx} + \lambda u = f	in \Omega = (a,b)
#	BC		u(a) = u(b) = 0 

# domain
(a,b)=(-1,1)
# subintervals and number of points per subinterval
S, Ns = 2, 15;
# 
a1, a2, a3 = -1,0,1

#x1,x2 = linspace(a1,a2,N), linspace(a2,a3,N)
l1,l2 = 1,1
#xi1, xi2 = 2./l1*(x1-a1)-1, 2./l2*(x2-a2)-1
xi1, xi2 = cos(pi/Ns*arange(0,Ns+1)), cos(pi/Ns*arange(0,Ns+1))
x1, x2 = l1/2*(xi1-1)+a2, l2/2*(xi2-1)+a3
xall = hstack((x1[:],x2[:]))
xm1 = hstack((x1[:-1],x2[1:]))
xm2 = hstack((x1[:-2],x2[2:]))

print(xall.shape)
print(xm1.shape)
print(xm2.shape)
plot(x1,0*x1,'o'), show()
