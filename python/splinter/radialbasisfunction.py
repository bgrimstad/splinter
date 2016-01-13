# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from .function import Function
from .utilities import *
from .datatable import DataTable


class RBFType:
    THIN_PLATE_SPLINE = 1
    MULTIQUADRIC = 2
    INVERSE_QUADRIC = 3
    INVERSE_MULTIQUADRIC = 4
    GAUSSIAN = 5


class RadialBasisFunction(Function):
    def __init__(self, dataOrFileName, rbfType=None, normalized=None):
        self._handle = None

        # If string we load the PolynomialRegression from the file
        if isString(dataOrFileName):
            fileName = getCString(dataOrFileName)
            self._handle = splinter._call(splinter._getHandle().rbf_load_init, fileName)

        # Else, construct PolynomialRegression from DataTable
        else:
            dataTable = DataTable(dataOrFileName)
            if rbfType is None:
                rbfType = RBFType.THIN_PLATE_SPLINE
            if rbfType < 1 or rbfType > 5:
                raise Exception("Error: Invalid rbfType: " + str(rbfType) + ". For valid values see the class RBFType.")

            if normalized is None:
                normalized = 0
            normalized = 0 if normalized == 0 else 1

            self._handle = splinter._call(splinter._getHandle().rbf_init, dataTable._getHandle(), rbfType, normalized)

        self._numVariables = splinter._call(splinter._getHandle().function_get_num_variables, self._handle)
