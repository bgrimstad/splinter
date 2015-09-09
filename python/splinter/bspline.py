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


class BSpline(Approximant):
	def __init__(self, dataTableOrFileName, degree=None):
		# If string we load the BSpline from the file
		if isString(dataTableOrFileName):
			fileName = getCString(dataTableOrFileName)
			self._handle = splinter._call(splinter._getHandle().bspline_load_init, fileName)
		
		# Else, construct BSpline from DataTable
		else:
			dataTable = dataTableOrFileName
			self._handle = splinter._call(splinter._getHandle().bspline_init, dataTable._getHandle(), degree)
			
		self._numVariables = splinter._call(splinter._getHandle().get_num_variables, self._handle)
