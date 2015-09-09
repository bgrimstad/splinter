# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import sys
from approximant import Approximant
from utilities import *


class RBFType:
	THIN_PLATE_SPLINE = 1
	MULTIQUADRIC = 2
	INVERSE_QUADRIC = 3
	INVERSE_MULTIQUADRIC = 4
	GAUSSIAN = 5

class RadialBasisFunction(Approximant):
	def __init__(self, dataTableOrFileName, rbfType=None):
		# If string we load the PolynomialRegression from the file
		if isString(dataTableOrFileName):
			fileName = getCString(dataTableOrFileName)
			self._handle = SPLINTER.call(SPLINTER.getHandle().rbf_load_init, fileName)
		
		# Else, construct PolynomialRegression from DataTable
		else:
			if rbfType is None:
				rbfType = RBFType.THIN_PLATE_SPLINE
			if rbfType < 1 or rbfType > 5:
				raise Exception("Error: Invalid rbfType: " + str(rbfType) + ". For valid values see the class RBFType.")
			
			dataTable = dataTableOrFileName
			self._handle = SPLINTER.call(SPLINTER.getHandle().rbf_init, dataTable.getHandle(), rbfType)
		
		
if __name__ == "__main__":
	import SPLINTER
	SPLINTER.load("/home/anders/SPLINTER/build/release/libsplinter-matlab-1-4.so")
	
	
	from datatable import DataTable
	
	def f(x):
		return x[0]*x[1]
	
	d = DataTable()
	for i in range(10):
		for j in range(10):
			d.addSample([i,j], f([i,j]))
	
	b = RadialBasisFunction(d, RBFType.GAUSSIAN)
	for i in range(10):
		for j in range(10):
			print(str(b.eval([0.9*i,0.9*j])) + " == " + str(0.81*i*j))
	
	print(b.evalJacobian([3,3]))
	print(b.evalHessian([3,3]))
	
	b.save("test.rbf")
	b2 = RadialBasisFunction("test.rbf")
	
	print(b.eval([2,3]))
	print(b2.eval([2,3]))