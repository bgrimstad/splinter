# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from ctypes import *
from .utilities import *


class Approximant:
	def __init__(self):
		self._handle = None
		self._numVariables = None

	def eval(self, x):
		if not isinstance(x, list):
			x = [x]
		
		# See if x is on the form [[x0,x1],[x2,x3]]
		# if not we assume it to be on the form
		# [x0, x1, x2, x3]
		if isinstance(x[0], list):
			x = flattenList(x)
		
		numPoints = len(x) // self._numVariables
		res = splinter._call(splinter._getHandle().eval, self._handle, (c_double * len(x))(*x), len(x))
		
		return CArrayToList(res, numPoints)

	def evalJacobian(self, x):
		if not isinstance(x, list):
			x = [x]
		
		# See if x is on the form [[x0,x1],[x2,x3]]
		# if not we assume it to be on the form
		# [x0, x1, x2, x3]
		if isinstance(x[0], list):
			x = flattenList(x)
		
		numPoints = len(x) // self._numVariables
		jac = splinter._call(splinter._getHandle().eval_jacobian, self._handle, (c_double * len(x))(*x), len(x))
		
		# Convert from ctypes array to Python list of lists
		# jacobians is a list of the jacobians in all evaluated points
		jacobians = []
		for i in range(numPoints):
			jacobians.append([])
			for j in range(self._numVariables):
				jacobians[i].append(jac[i*self._numVariables + j])
		return jacobians

	def evalHessian(self, x):
		if not isinstance(x, list):
			x = [x]
		
		# See if x is on the form [[x0,x1],[x2,x3]]
		# if not we assume it to be on the form
		# [x0, x1, x2, x3]
		if isinstance(x[0], list):
			x = flattenList(x)
		
		numPoints = len(x) // self._numVariables
		hes = splinter._call(splinter._getHandle().eval_hessian, self._handle, (c_double * len(x))(*x), len(x))
		
		# Convert from ctypes array to Python list of list of lists
		# hessians is a list of the hessians in all points
		hessians = []
		for i in range(numPoints):
			hessians.append([])
			for j in range(self._numVariables):
				hessians[i].append([])
				for k in range(self._numVariables):
					hessians[i][j].append(hes[i*self._numVariables*self._numVariables + j*self._numVariables + k])
		return hessians

	def getNumVariables(self):
		return splinter._call(splinter._getHandle().get_num_variables, self._handle)

	def save(self, fileName):
		splinter._call(splinter._getHandle().save, self._handle, getCString(fileName))

	def __del__(self):
		if self._handle is not None:
			splinter._call(splinter._getHandle().delete_approximant, self._handle)
		self._handle = None
