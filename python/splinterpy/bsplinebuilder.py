# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from .splinter_backend import splinter_backend_obj as backend
from .bspline import BSpline
from .datatable import DataTable
from typing import Union, List
from .utilities import *


Smoothing = BSpline.Smoothing
KnotSpacing = BSpline.KnotSpacing


class BSplineBuilder:

    def __init__(self, dim_x, dim_y, degree: int=3, knot_spacing: int=KnotSpacing.AS_SAMPLED,
                 num_basis_functions: int=int(1e6)):
        self._handle = None  # Handle for referencing the c side of this object
        self._dim_x = dim_x
        self._dim_y = dim_y
        self._num_basis_functions = [10 ** 3] * self._dim_x

        self._degrees = None
        self._knot_spacing = None
        self._num_basis_functions = None

        f_handle_init = backend.handle.splinter_bspline_builder_init
        self._handle = backend.call(f_handle_init, self._dim_x, self._dim_y)
        self.degree(degree)
        self.knot_spacing(knot_spacing)
        self.num_basis_functions(num_basis_functions)

    def degree(self, degrees: Union[List[int], int]) -> 'BSplineBuilder':
        # If the value is a single number, make it a list of numVariables length
        if not isinstance(degrees, list):
            degrees = [degrees] * self._dim_x

        if len(degrees) != self._dim_x:
            raise ValueError("Inconsistent number of degrees.")

        valid_degrees = range(0, 6)
        for deg in degrees:
            if deg not in valid_degrees:
                raise ValueError("Invalid degree: " + str(deg))

        self._degrees = degrees

        f_handle = backend.handle.splinter_bspline_builder_set_degree
        backend.call(f_handle, self._handle, list_to_c_array_of_ints(self._degrees), len(self._degrees))
        return self

    def knot_spacing(self, knot_spacing: int) -> 'BSplineBuilder':
        if not KnotSpacing.is_valid(knot_spacing):
            raise ValueError("Invalid knotspacing: " + str(knot_spacing))

        self._knot_spacing = knot_spacing

        f_handle = backend.handle.splinter_bspline_builder_set_knot_spacing
        backend.call(f_handle, self._handle, self._knot_spacing)
        return self

    def num_basis_functions(self, num_basis_functions: Union[List[int], int]) -> 'BSplineBuilder':
        # If the value is a single number, make it a list of num_variables length
        if not isinstance(num_basis_functions, list):
            num_basis_functions = [num_basis_functions] * self._dim_x

        if len(num_basis_functions) != self._dim_x:
            raise ValueError("Inconsistent number of degrees.")

        for num_basis_func in num_basis_functions:
            if not isinstance(num_basis_func, int):
                raise TypeError("Number of basis functions not integer: " + str(
                        num_basis_func))
            if num_basis_func < 1:
                raise ValueError("Number of basis functions < 1: " + str(
                        num_basis_func))

        self._num_basis_functions = num_basis_functions

        f_handle = backend.handle.splinter_bspline_builder_set_num_basis_functions
        backend.call(f_handle, self._handle, list_to_c_array_of_ints(self._num_basis_functions),
                     len(self._num_basis_functions))
        return self

    # Returns a handle to the created internal BSpline object
    def fit(self, X, Y, smoothing: int=Smoothing.NONE, alpha: float=0.1, weights: List[float]=list()) -> BSpline:

        if alpha < 0:
            raise ValueError("'alpha' must be non-negative")

        if not Smoothing.is_valid(smoothing):
            raise ValueError("Invalid smoothing type: " + str(smoothing))

        data_table = DataTable(X, Y)
        f_handle = backend.handle.splinter_bspline_builder_fit
        bspline_handle = backend.call(f_handle, self._handle, data_table._get_handle(), smoothing, alpha,
                                      list_to_c_array_of_doubles(weights), len(weights))
        return BSpline(bspline_handle)

    def __del__(self):
        if self._handle is not None:
            f_handle = backend.handle.splinter_bspline_builder_delete
            backend.call(f_handle, self._handle)
        self._handle = None


def bspline_interpolator(X, Y, degree: int=3) -> BSpline:

    if degree < 0:
        raise ValueError("'degree' must be non-negative")

    data_table = DataTable(X, Y)
    f_handle = backend.handle.splinter_bspline_interpolator
    bspline_handle = backend.call(f_handle, data_table._get_handle(), degree)
    return BSpline(bspline_handle)


def bspline_smoother(X, Y, degree: int=3, smoothing: Smoothing=Smoothing.PSPLINE, alpha: float=0.1,
                     weights: List[float]=list()) -> BSpline:

    if degree < 0:
        raise ValueError("'degree' must be non-negative")

    if alpha < 0:
        raise ValueError("'alpha' must be non-negative")

    data_table = DataTable(X, Y)
    f_handle = backend.handle.splinter_bspline_smoother
    bspline_handle = backend.call(f_handle, data_table._get_handle(), degree, smoothing, alpha,
                                  list_to_c_array_of_doubles(weights), len(weights))
    return BSpline(bspline_handle)


def bspline_unfitted(X, Y, degrees: List[int], knot_spacing: KnotSpacing, num_basis_functions: List[int]) -> BSpline:

    if not isinstance(degrees, list):
        raise ValueError("'degrees' must be of type list")

    for d in degrees:
        if d < 0:
            raise ValueError("Degrees must be non-negative")

    if not KnotSpacing.is_valid(knot_spacing):
        raise ValueError("Invalid knot spacing type (see KnotSpacing)")

    if not isinstance(num_basis_functions, list):
        raise ValueError("'num_basis_functions' must be of type list")

    if len(degrees) != len(num_basis_functions):
        raise ValueError("Length of 'degrees' and 'num_basis_functions' must be equal")

    for d, num in zip(degrees, num_basis_functions):
        if num < d + 1:
            raise ValueError("Number of basis functions must be at least degree + 1")

    f_handle = backend.handle.splinter_bspline_unfitted
    data_table = DataTable(X, Y)
    _degrees = list_to_c_array_of_ints(degrees)
    _num_basis_functions = list_to_c_array_of_ints(num_basis_functions)
    bspline_handle = backend.call(f_handle, data_table._get_handle(), _degrees, len(degrees), knot_spacing,
                                  _num_basis_functions, len(num_basis_functions))
    return BSpline(bspline_handle)

