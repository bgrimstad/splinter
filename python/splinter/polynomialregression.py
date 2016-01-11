# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from .function import Function
from .utilities import *


class PolynomialRegression(Function):
    def __init__(self, dataTableOrFileName, degree=None):
        # If string we load the PolynomialRegression from the file
        if isString(dataTableOrFileName):
            fileName = getCString(dataTableOrFileName)
            self._handle = splinter._call(splinter._getHandle().polynomial_regression_load_init, fileName)

        # Else, construct PolynomialRegression from DataTable
        else:
            if not isinstance(degree, list):
                degree = [degree]*dataTableOrFileName.getNumVariables()
            dataTable = dataTableOrFileName
            self._handle = splinter._call(splinter._getHandle().polynomial_regression_init, dataTable._getHandle(), listToCArrayOfInts(degree), len(degree))

        self._numVariables = splinter._call(splinter._getHandle().function_get_num_variables, self._handle)
