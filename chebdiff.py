#!/usr/bin/python
from numpy import *
from pylab import plot,show,figure,savefig,semilogy

def cheb(N,a,b):
	y=-cos(linspace(0,pi,N+1))
	x=.5*(a+b + (b-a)*y)
	D=zeros((N+1,N+1))
	for ii in range(N+1):
		for jj in range(N+1):
			if ii==jj:
				D[jj,jj]=sum([1./(x[jj]-x[k]) for k in range(N+1) if k != jj])
			else:
				ai = prod([(x[ii] - x[k]) for k in range(N+1) if k != ii])
				aj = prod([(x[jj] - x[k]) for k in range(N+1) if k != jj])
				D[ii,jj]=ai/aj*1./(x[ii]-x[jj])
	return (D,x)
	
if __name__ == '__main__':
	N=30;
	(D,x) = cheb(N,-1,1)
	u=exp(x)*sin(5*x)
	ux=dot(D,u)
	ux_exact=exp(x)*(sin(5*x)+5*cos(5*x))
	error = ux - ux_exact
	# test domain decomposition
	(D1,x1)=cheb(N/2,-1,0)
	u1=exp(x1)*sin(5*x1)
	u1x=dot(D1,u1)
	u1x_exact=exp(x1)*(sin(5*x1)+5*cos(5*x1))
	(D2,x2)=cheb(N/2,0,1)
	u2=exp(x2)*sin(5*x2)
	u2x=dot(D2,u2)
	u2x_exact=exp(x2)*(sin(5*x2)+5*cos(5*x2))
	figure(1);plot(x1,u1x,'.-',x1,u1x_exact,x2,u2x,'.-',x2,u2x_exact)
	savefig('comparison_domains.png')
	figure(2);semilogy(x1,abs(u1x-u1x_exact),x2,abs(u2x-u2x_exact))
	savefig('comparison_domains_err.png')
	figure(3);plot(hstack((x1,x2)),0*hstack((x1,x2)),'o')
	savefig('cheb_points.png');