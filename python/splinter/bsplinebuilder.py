# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from ctypes import *
from .bspline import BSpline
from .datatable import DataTable
from typing import Union, List


class BSplineBuilder:

    class Smoothing:
        NONE, IDENTITY, PSPLINE = range(3)

        @staticmethod
        def is_valid(value):
            return value in range(3)

    class KnotSpacing:
        AS_SAMPLED, EQUIDISTANT, EXPERIMENTAL = range(3)

        @staticmethod
        def is_valid(value):
            return value in range(3)

    def __init__(self, x, y, degree: int=3, smoothing: int=Smoothing.NONE, alpha: float=0.1,
                 knot_spacing: int=KnotSpacing.AS_SAMPLED, num_basis_functions: int=int(1e6)):
        self._handle = None  # Handle for referencing the c side of this object
        self._datatable = DataTable(x, y)
        self._num_basis_functions = [10 ** 3] * self._datatable.get_num_variables()

        self._degrees = None
        self._alpha = None
        self._smoothing = None
        self._knot_spacing = None
        self._num_basis_functions = None

        f_handle_init = splinter._get_handle().splinter_bspline_builder_init
        self._handle = splinter._call(f_handle_init, self._datatable._get_handle())
        self.degree(degree)
        self.set_alpha(alpha)
        self.smoothing(smoothing)
        self.knot_spacing(knot_spacing)
        self.num_basis_functions(num_basis_functions)

    def degree(self, degrees: Union[List[int], int]) -> 'BSplineBuilder':
        # If the value is a single number, make it a list of numVariables length
        if not isinstance(degrees, list):
            degrees = [degrees] * self._datatable.get_num_variables()

        if len(degrees) != self._datatable.get_num_variables():
            raise ValueError("Inconsistent number of degrees.")

        valid_degrees = range(0, 6)
        for deg in degrees:
            if deg not in valid_degrees:
                raise ValueError("Invalid degree: " + str(deg))

        self._degrees = degrees

        f_handle = splinter._get_handle().splinter_bspline_builder_set_degree
        splinter._call(f_handle, self._handle, (c_int * len(self._degrees))(*self._degrees), len(self._degrees))
        return self

    def set_alpha(self, new_alpha: float) -> 'BSplineBuilder':
        if new_alpha < 0:
            raise ValueError("'alpha' must be non-negative.")

        self._alpha = new_alpha

        f_handle = splinter._get_handle().splinter_bspline_builder_set_alpha
        splinter._call(f_handle, self._handle, self._alpha)
        return self

    def smoothing(self, smoothing: int) -> 'BSplineBuilder':
        if not BSplineBuilder.Smoothing.is_valid(smoothing):
            raise ValueError("Invalid smoothing: " + str(smoothing))

        self._smoothing = smoothing

        f_handle = splinter._get_handle().splinter_bspline_builder_set_smoothing
        splinter._call(f_handle, self._handle, self._smoothing)
        return self

    def knot_spacing(self, knot_spacing: int) -> 'BSplineBuilder':
        if not BSplineBuilder.KnotSpacing.is_valid(knot_spacing):
            raise ValueError("Invalid knotspacing: " + str(knot_spacing))

        self._knot_spacing = knot_spacing

        f_handle = splinter._get_handle().splinter_bspline_builder_set_knot_spacing
        splinter._call(f_handle, self._handle, self._knot_spacing)
        return self

    def num_basis_functions(self, num_basis_functions: Union[List[int], int]) -> 'BSplineBuilder':
        # If the value is a single number, make it a list of num_variables length
        if not isinstance(num_basis_functions, list):
            num_basis_functions = [num_basis_functions] * self._datatable.get_num_variables()

        if len(num_basis_functions) != self._datatable.get_num_variables():
            raise ValueError("Inconsistent number of degrees.")

        for num_basis_func in num_basis_functions:
            if not isinstance(num_basis_func, int):
                raise TypeError("Number of basis functions not integer: " + str(
                        num_basis_func))
            if num_basis_func < 1:
                raise ValueError("Number of basis functions < 1: " + str(
                        num_basis_func))

        self._num_basis_functions = num_basis_functions

        f_handle = splinter._get_handle().splinter_bspline_builder_set_num_basis_functions
        splinter._call(f_handle, self._handle, (c_int * len(self._num_basis_functions))(*self._num_basis_functions),
                       len(self._num_basis_functions))
        return self

    # Returns a handle to the created internal BSpline object
    def build(self) -> BSpline:
        f_handle = splinter._get_handle().splinter_bspline_builder_build
        bspline_handle = splinter._call(f_handle, self._handle)

        return BSpline(bspline_handle)

    def __del__(self):
        if self._handle is not None:
            f_handle = splinter._get_handle().splinter_bspline_builder_delete
            splinter._call(f_handle, self._handle)
        self._handle = None
