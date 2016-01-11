# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from .datatable import DataTable
from ctypes import *
from .bspline import BSpline


class BSplineBuilder:
    class Degree:
        LINEAR, QUADRATIC, CUBIC, QUARTIC = range(1, 5)

        @staticmethod
        def isValid(value):
            return value in range(1, 5)

    class Smoothing:
        NONE, REGULARIZATION, PSPLINE = range(3)

        @staticmethod
        def isValid(value):
            return value in range(3)

    class KnotSpacing:
        SAMPLE, EQUIDISTANT, EXPERIMENTAL = range(3)

        @staticmethod
        def isValid(value):
            return value in range(3)

    def __init__(self, datatable):
        self._handle = None  # Handle for referencing the c side of this object

        self._datatable = datatable
        self._degrees = [BSplineBuilder.Degree.CUBIC] * datatable.getNumVariables()
        self._numBasisFunctions = [10**3] * datatable.getNumVariables()
        self._knotSpacing = BSplineBuilder.KnotSpacing.SAMPLE
        self._smoothing = BSplineBuilder.Smoothing.NONE
        self._lambda = 0.03

        self._handle = splinter._call(splinter._getHandle().bsplinebuilder_init, datatable._getHandle())

    def degree(self, degrees):
        # If the value is a single number, make it a list of numVariables length
        if not isinstance(degrees, list):
            degrees = [degrees] * self._datatable.getNumVariables()

        if len(degrees) != self._datatable.getNumVariables():
            raise ValueError("BSplineBuilder:degree: Inconsistent number of degrees.")

        for degree in degrees:
            if not BSplineBuilder.Degree.isValid(degree):
                raise ValueError("BSplineBuilder:degree: Invalid degree: " + str(degree))

        self._degrees = degrees

        splinter._call(splinter._getHandle().bsplinebuilder_set_degree, self._handle, (c_int * len(self._degrees))(*self._degrees), len(self._degrees))
        return self

    def numBasisFunctions(self, numBasisFunctions):
        # If the value is a single number, make it a list of numVariables length
        if not isinstance(numBasisFunctions, list):
            numBasisFunctions = [numBasisFunctions] * self._datatable.getNumVariables()

        if len(numBasisFunctions) != self._datatable.getNumVariables():
            raise ValueError("BSplineBuilder:numBasisFunctions: Inconsistent number of degrees.")

        for numBasisFunction in numBasisFunctions:
            if not isinstance(numBasisFunction, int):
                raise ValueError("BSplineBuilder:numBasisFunctions: Invalid number of basis functions (must be integer): " + str(numBasisFunction))

        self._numBasisFunctions = numBasisFunctions

        splinter._call(splinter._getHandle().bsplinebuilder_set_num_basis_functions, self._handle, (c_int * len(self._numBasisFunctions))(*self._numBasisFunctions), len(self._numBasisFunctions))
        return self

    def knotSpacing(self, knotSpacing):
        if not BSplineBuilder.KnotSpacing.isValid(knotSpacing):
            raise ValueError("BSplineBuilder::knotSpacing: Invalid knotspacing: " + str(knotSpacing))

        self._knotSpacing = knotSpacing

        splinter._call(splinter._getHandle().bsplinebuilder_set_knot_spacing, self._handle, self._knotSpacing)
        return self

    def smoothing(self, smoothing):
        if not BSplineBuilder.KnotSpacing.isValid(smoothing):
            raise ValueError("BSplineBuilder::smoothing: Invalid smoothing: " + str(smoothing))

        self._smoothing = smoothing

        splinter._call(splinter._getHandle().bsplinebuilder_set_smoothing, self._handle, self._smoothing)
        return self

    def setLambda(self, newLambda):
        if newLambda < 0:
            raise ValueError("BSplineBuilder:setLambda: Lambda must be non-negative.")

        self._lambda = newLambda

        splinter._call(splinter._getHandle().bsplinebuilder_set_lambda, self._handle, self._lambda)
        return self

    # Returns a handle to the created internal BSpline object
    def build(self):
        bspline_handle = splinter._call(splinter._getHandle().bsplinebuilder_build, self._handle)

        return BSpline(bspline_handle)
