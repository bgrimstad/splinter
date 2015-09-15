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


class PSpline(Approximant):
	def __init__(self, dataTableOrFileName, smoothing = None):
		# If string we load the PSpline from the file
		if isString(dataTableOrFileName):
			fileName = getCString(dataTableOrFileName)
			self._handle = splinter._call(splinter._getHandle().pspline_load_init, fileName)
		
		# Else, construct PSpline from DataTable
		else:
			dataTable = dataTableOrFileName
			if smoothing is None:
				smoothing = 0.03
			self._handle = splinter._call(splinter._getHandle().pspline_init, dataTable._getHandle(), smoothing)
			
		self._numVariables = splinter._call(splinter._getHandle().approximant_get_num_variables, self._handle)
