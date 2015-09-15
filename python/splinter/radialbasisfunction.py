# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import sys
from . import splinter
from .approximant import Approximant
from .utilities import *


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
			self._handle = splinter._call(splinter._getHandle().rbf_load_init, fileName)
		
		# Else, construct PolynomialRegression from DataTable
		else:
			if rbfType is None:
				rbfType = RBFType.THIN_PLATE_SPLINE
			if rbfType < 1 or rbfType > 5:
				raise Exception("Error: Invalid rbfType: " + str(rbfType) + ". For valid values see the class RBFType.")
			
			dataTable = dataTableOrFileName
			self._handle = splinter._call(splinter._getHandle().rbf_init, dataTable._getHandle(), rbfType)
			
		self._numVariables = splinter._call(splinter._getHandle().get_num_variables, self._handle)
