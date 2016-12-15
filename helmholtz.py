#!/usr/bin/python
from __future__ import division
from pylab import *
from numpy import *
from time import sleep

# solve 
#	PDE		-u_{xx} + \lambda u = f	in \Omega = (a,b)
#	BC		u(a) = u(b) = 0 

# domain
(a,b)=(-1,1)
# subintervals and number of points per subinterval
S, Ns = 2, 14;
# 
a1, a2, a3 = -1,0,1

#x1,x2 = linspace(a1,a2,N), linspace(a2,a3,N)
l1,l2 = 1,1
#xi1, xi2 = 2./l1*(x1-a1)-1, 2./l2*(x2-a2)-1
xi1, xi2 = cos(pi/Ns*arange(0,Ns+1)), cos(pi/Ns*arange(0,Ns+1))
x1, x2 = l1/2*(xi1+1)+a1, l2/2*(xi2+1)+a2
x = hstack((x1[1:],x2[:-1]))
xi = hstack((xi1,xi2))
dx=x[1]-x[0]
# T_n = cos(n*arccos(x))
T=lambda n,x: cos(n*arccos(x))
cbar=lambda n:  1 + (n==0 or n==Ns)*1

def J(k):
	if k==0:
		return 0
	else:
		return -4*sum( [1/(2*q-1) for q in range(1,k+1)] )
		
def amat(n,m):
	if (n+m) % 2==1:
		return 0
	else:
		return n*m/2*(J(abs(n-m)//2) - J(abs(n+m)//2))
		
def bmat(n,m):
	if (n+m) % 2==1:
		return 0
	else:
		return 1/(1-(n+m)**2) + 1/(1-(n-m)**2)
		
def Amat(j,k,s):
	psum = 0
	for n in range(Ns+1):
		for m in range(Ns+1):
			psum += 1/(cbar(n)*cbar(m))*T(n,xi[k])*T(m,xi[j])*amat(n,m)
	return 8/((Ns)**2*cbar(k)*cbar(j))*psum
	
def Bmat(j,k,s):
	psum = 0
	for n in range(Ns+1):
		for m in range(Ns+1):
			psum += 1/(cbar(n)*cbar(m))*T(n,xi[k])*T(m,xi[j])*bmat(n,m)
	return 2*l1/((Ns**2)*cbar(k)*cbar(j))*psum

lamb = 0;
dim = (S*Ns,S*Ns);
A,B,C = zeros(dim), zeros(dim), zeros(dim)
for i in range(A.shape[0]):
	for j in range(A.shape[1]):
		if i//Ns == j//Ns:		# ha! clever!
			A[i,j] = Amat(i,j,(i)//Ns)
			B[i,j] = Bmat(i,j,(i)//Ns)
# condense A matrix
#A[x==-1,:],A[x==1,:],A[:,x==-1],A[:,x==1]=0,0,0,0
# But then just set A[x==\pm 1]=1 so that we don't have to chop matrix
A[x==-1,x==-1],A[x==1, x==1 ] = 1, 1
# now pair x=dx with x=-dx
#A[0,:], A[-1,:] = 0,0
A[0,-1],A[-1,0]=Amat(0,-1,1),Amat(-1,0,1)
#for rowelt in range(A.shape[1]):
#	A[0,rowelt] = Amat(0,rowelt,1)
#for rowelt in range(A.shape[1]):
#	A[A.shape[0]-1,rowelt]=Amat(A.shape[0],rowelt,0)
# ditto, for B matrix
B[x==-1,:],B[x==1,:],B[:,x==-1],B[:,x==1]=0,0,0,0
B[x==-1,x==-1],B[x==1, x==1 ] = 1, 1
#B[0,-1],B[-1,0]=1,1


C = A + lamb * B
f = cos(pi*x + pi/4)
f[x==1]=0
f[x==-1]=0
uexact = (sqrt(2)+sqrt(2)*cos(pi*x)-sqrt(2)*sin(pi*x))/(2*pi**2)
uexactm = -(sqrt(2)+sqrt(2)*cos(pi*x)-sqrt(2)*sin(pi*x))/(2*pi**2)

figure(1),spy(A),draw()
#exit()
rhs = dot(B,f)
u = linalg.solve(C,rhs)
figure(2),plot(x,f,'o');legend(("Starting function",)), draw()
figure(3),plot(x,u,'o',x,uexact,'.',x,uexactm,'.');legend(("Second Deriv","Exact","Exact * (-1)"), loc='lower right'); show()