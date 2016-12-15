#!/usr/bin/python
# dgsem1d.py
#
# Created by Travis Johnson on 2010-06-01.
# Copyright (c) 2010 . All rights reserved.
from __future__ import division
from pylab import *
from numpy import *
from KoprivaMethods import *

class Element:
	def __init__(self, dG, nEqn, xL, xR):
		"""docstring for __init__"""
		self.dx = 
		self.xL = xL
		self.xR = xR
		self.nEqn = nEqn
		self.Q = zeros((N,nEqn))
		self.Qdot = zeros((N,nEqn))
		self.G = zeros((N,nEqn))
		self.QL = zeros((nEqn, 1))
		self.QR = zeros((nEqn, 1))
		self.FStarL = zeros((nEqn,1))
		self.FStarR = zeros((nEqn,1))
	def InterpolateToBoundaries(self):
		"""docstring for InterpolateToBoundaries"""
		pass
	def LocalTimeDerivative():
		"""docstring for LocalTimeDerivative"""
		pass
	def AffineMap(self, xi):
		"""docstring for AffineMap"""
		pass
	
	

class Mesh1D:
	def __init__(self, K, N, x):
		"""docstring for __init__"""
		self.K = K
		self.eltk = []
		self.pk = []
	def GlobalTimeDerivative(self):
		"""docstring for GlobalTimeDerivative"""
		pass
	


