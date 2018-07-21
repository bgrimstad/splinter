# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from .splinter_backend import splinter_backend_obj as backend
from .bspline import BSpline
from .datatable import DataTable
from typing import List
from .utilities import *


Smoothing = BSpline.Smoothing
KnotSpacing = BSpline.KnotSpacing


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

