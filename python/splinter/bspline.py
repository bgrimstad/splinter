# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from .function import Function
from .utilities import *


class BSpline(Function):
    def __init__(self, handleOrFileName):
        super(BSpline, self).__init__()

        # If string we load the BSpline from the file
        if isString(handleOrFileName):
            fileName = getCString(handleOrFileName)
            self._handle = splinter._call(splinter._getHandle().splinter_bspline_load_init, fileName)

        # Else, the argument is the handle to the internal BSpline object
        else:
            self._handle = handleOrFileName

        self._numVariables = splinter._call(splinter._getHandle().splinter_bspline_get_num_variables, self._handle)

    def getKnotVectors(self):
        """
        :return List of knot vectors (of possibly differing lengths)
        """
        # Get the sizes of the knot vectors first
        knotVectorSizesRaw = splinter._call(splinter._getHandle().splinter_bspline_get_knot_vector_sizes, self._handle)
        knotVectorSizes = CArrayToList(knotVectorSizesRaw, self._numVariables)

        knotVectorsRaw = splinter._call(splinter._getHandle().splinter_bspline_get_knot_vectors, self._handle)
        # Get the knot vectors as one long vector where knot vectors v1, ..., vn is laid out like this:
        # v11, ..., v1m, ..., vn1, ..., vno
        totSize = sum(knotVectorSizes)
        knotVectorsSerialized = CArrayToList(knotVectorsRaw, totSize)

        # Then reconstructor the knot vectors from the long vector by utilizing that we know how long each vector is
        knotVectors = []
        start = 0
        for knotVectorSize in knotVectorSizes:
            knotVectors.append(knotVectorsSerialized[start:start+knotVectorSize])
            start += knotVectorSize

        return knotVectors

    def getCoefficients(self):
        """
        :return List of the coefficients of the BSpline
        """
        numCoefficients = splinter._call(splinter._getHandle().splinter_bspline_get_num_coefficients, self._handle)
        coefficientsRaw = splinter._call(splinter._getHandle().splinter_bspline_get_coefficients, self._handle)

        return CArrayToList(coefficientsRaw, numCoefficients)

    def getControlPoints(self):
        """
        Get the matrix with the control points of the BSpline.
        :return Matrix (as a list of lists) with getNumVariables+1 columns and len(getCoefficients) rows
        """
        controlPointsRaw = splinter._call(splinter._getHandle().splinter_bspline_get_control_points, self._handle)

        # Yes, num_coefficients is correct
        numRows = splinter._call(splinter._getHandle().splinter_bspline_get_num_coefficients, self._handle)
        numCols = self._numVariables+1

        controlPointsFlattened = CArrayToList(controlPointsRaw, numRows*numCols)

        controlPoints = []
        start = 0
        for row in range(numRows):
            controlPoints.append(controlPointsFlattened[start:start+numCols])
            start += numCols

        return controlPoints
