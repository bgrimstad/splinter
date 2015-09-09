# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import SPLINTER
from ctypes import *
from utilities import *


class Approximant:
	def __init__(self):
		self._handle = None
	
	def eval(self, x):
		return SPLINTER.call(SPLINTER.getHandle().eval, self._handle, (c_double * len(x))(*x), len(x))
	
	def evalJacobian(self, x):
		jac = SPLINTER.call(SPLINTER.getHandle().eval_jacobian, self._handle, (c_double * len(x))(*x), len(x))
		# Convert from ctypes array to Python list
		return [jac[i] for i in range(len(x))]
	
	def evalHessian(self, x):
		hes = SPLINTER.call(SPLINTER.getHandle().eval_hessian, self._handle, (c_double * len(x))(*x), len(x))
		# Convert from ctypes array to Python list of lists
		hesList = []
		for i in range(len(x)):
			hesList.append([])
			for j in range(len(x)):
				hesList[i].append(hes[i*len(x) + j])
		
		return hesList
	
	def getNumVariables(self):
		return SPLINTER.call(SPLINTER.getHandle().get_num_variables, self._handle)
	
	def save(self, fileName):
		SPLINTER.call(SPLINTER.getHandle().save, self._handle, getCString(fileName))
		
	def load(self, fileName):
		SPLINTER.call(SPLINTER.getHandle().load, self._handle, getCString(fileName))
		
	def __del__(self):
		if self._handle is not None:
			SPLINTER.call(SPLINTER.getHandle().delete_approximant, self._handle)
		self._handle = None