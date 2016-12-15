#!/usr/bin/python
# deriv_errs.py
#
# Created by Travis Johnson on 2010-06-02.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *

f = lambda x: exp(2*x)*sin(6*x)
exact = lambda x: 2*exp(2*x)*sin(6*x) + 6*exp(2*x)*cos(6*x)
f = lambda x: 3/(5-4*cos(x))
exact = lambda x: -12*sin(x)/(5-4*cos(x))**2
#f = lambda x: 
nVec = range(6,36,2)
errLG, errLGL, errCh, errCGL = [],[],[],[]
maxeigLG, maxeigLGL, maxeigCh, maxeigCGL = [],[],[],[]
for N in nVec:
	xLG, wLG = LegendreGaussNodesAndWeights(N)
	xLGL,wLGL = LegendreGaussLobattoNodesAndWeights(N)
	xCh, wCh = ChebyshevGaussNodesAndWeights(N)
	xCGL,wCGL = ChebyshevGaussLobattoNodesAndWeights(N)

	dLG = PolynomialDerivativeMatrix(xLG)
	#figure(1),plot(xLG, exact(xLG),xLG, dot(dLG,f(xLG)),'o')
	errLG.append(linalg.norm(exact(xLG)-dot(dLG,f(xLG))))
	maxeigLG.append(max(abs(eigvals(dLG[1:,1:]))))
	
	dLGL = PolynomialDerivativeMatrix(xLGL)
	#figure(2),plot(xLGL, exact(xLGL),xLGL, dot(dLGL,f(xLGL)),'o')
	errLGL.append(linalg.norm(exact(xLGL)-dot(dLGL,f(xLGL))))
	maxeigLGL.append(max(abs(eigvals(dLGL[1:,1:]))))
	
	dCh = PolynomialDerivativeMatrix(xCh)
	#figure(3),plot(xCh, exact(xCh),xCh, dot(dCh,f(xCh)),'o')
	errCh.append(linalg.norm(exact(xCh)-dot(dCh,f(xCh))))
	maxeigCh.append(max(abs(eigvals(dCh[1:,1:]))))
	
	dCGL = PolynomialDerivativeMatrix(xCGL)
	#figure(4),plot(xCGL, exact(xCGL),xCGL, dot(dCGL,f(xCGL)),'o')
	errCGL.append(linalg.norm(exact(xCGL)-dot(dCGL,f(xCGL))))
	maxeigCGL.append(max(abs(eigvals(dCGL[1:,1:]))))
	
figure(1),semilogy(nVec, errLG, nVec, errLGL, nVec, errCh, nVec, errCGL)
title('Semi-Log of Error for various methods at several orders')
legend(('LG','LGL','CH','CGL'))
savefig('errors.eps')

print maxeigLG
print maxeigLGL
print maxeigCh
print maxeigCGL
figure(2),loglog(nVec, maxeigLG, nVec, maxeigLGL, nVec, maxeigCh, nVec, maxeigCGL)
title('Log-Log of eigenvalues for various methods at several orders')
legend(('LG','LGL','CH','CGL'),loc='lower right')
axis([min(nVec),max(nVec),min(maxeigLG),max(maxeigLG)])
savefig('eigengrowth.eps')


show()