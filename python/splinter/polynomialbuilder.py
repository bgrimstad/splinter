# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from ctypes import *
from .polynomial import Polynomial
from .builderbase import BuilderBase
from .utilities import flattenList


class PolynomialBuilder(BuilderBase):
    def __init__(self, data):
        super(PolynomialBuilder, self).__init__(data)

        self._degrees = [0] * self._datatable.getNumVariables()

        self._handle = splinter._call(splinter._getHandle().splinter_polynomial_builder_init, self._datatable._getHandle())

    def powers(self, powers):
        self._powers = flattenList(powers)

        splinter._call(splinter._getHandle().splinter_polynomial_builder_set_powers, self._handle, (c_int * len(self._powers))(*self._powers), len(self._powers))
        return self

    def setLambda(self, newLambda):
        if newLambda < 0:
            raise ValueError("PolynomialBuilder:lambda: Lambda must be non-negative.")

        self._lambda = newLambda

        splinter._call(splinter._getHandle().splinter_polynomial_builder_set_lambda, self._handle, self._lambda)
        return self

    # Returns a handle to the created internal Polynomial object
    def build(self):
        polynomial_handle = splinter._call(splinter._getHandle().splinter_polynomial_builder_build, self._handle)

        return Polynomial(polynomial_handle)

    def __del__(self):
        if self._handle is not None:
            splinter._call(splinter._getHandle().splinter_polynomial_builder_delete, self._handle)
        self._handle = None
