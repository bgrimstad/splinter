# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from .function import Function
from .utilities import *
from typing import List, Union

ListList = List[List[float]]
ControlPointsType = Union[ListList, List[float]]


class BSpline(Function):
    def __init__(self, handle_or_filename):
        super().__init__()

        # If string we load the BSpline from the file
        if is_string(handle_or_filename):
            filename = get_c_string(handle_or_filename)
            self._handle = splinter._call(splinter._get_handle().splinter_bspline_load_init, filename)

        # Else, the argument is the handle to the internal BSpline object
        else:
            self._handle = handle_or_filename

        self._dim_x = splinter._call(splinter._get_handle().splinter_bspline_get_dim_x, self._handle)
        self._dim_y = splinter._call(splinter._get_handle().splinter_bspline_get_dim_y, self._handle)

    @staticmethod
    def init_from_param(control_points: ControlPointsType, knot_vectors: ListList, degrees: List[int]) -> 'BSpline':
        """
        Builds B-splines given explicit parameters
        :param control_points: a list of control points in R^dim_y
        :param knot_vectors: a list of dim_x knot vectors of variable size
        :param degrees: a list of dim_x degrees (each being a non-negative integer)
        :return: BSpline object
        """
        if len(knot_vectors) != len(degrees):
            raise ValueError("Inconsistent data: len(knot_vectors) should equal len(degrees)")

        dim_x = len(degrees)

        if not any(control_points):
            raise ValueError("Array of control points is empty")

        # Handle 1-D case - create a list of lists from a list of control points
        if not any(isinstance(cp, list) for cp in control_points):
            control_points = [[cp] for cp in control_points]

        dim_y = len(control_points[0])

        control_points_c_array = list_to_c_array_of_doubles(flatten_list(control_points))
        num_control_points = len(control_points)
        knot_vectors_c_array = list_to_c_array_of_doubles(flatten_list(knot_vectors))
        num_knots_per_vector_c_array = list_to_c_array_of_ints(list([len(vec) for vec in knot_vectors]))
        degrees_c_array = list_to_c_array_of_ints(degrees)
        handle = splinter._call(splinter._get_handle().splinter_bspline_param_init,
                                dim_x,
                                dim_y,
                                control_points_c_array,
                                num_control_points,
                                knot_vectors_c_array,
                                num_knots_per_vector_c_array,
                                degrees_c_array)
        return BSpline(handle)

    def get_knot_vectors(self) -> ListList:
        """
        :return List of knot vectors (of possibly differing lengths)
        """
        # Get the sizes of the knot vectors first
        knot_vector_sizes_raw = splinter._call(splinter._get_handle().splinter_bspline_get_knot_vector_sizes, self._handle)
        knot_vector_sizes = c_array_to_list(knot_vector_sizes_raw, self._num_variables)

        knot_vectors_raw = splinter._call(splinter._get_handle().splinter_bspline_get_knot_vectors, self._handle)
        # Get the knot vectors as one long vector where knot vectors v1, ..., vn is laid out like this:
        # v11, ..., v1m, ..., vn1, ..., vno
        tot_size = sum(knot_vector_sizes)
        knot_vectors_serialized = c_array_to_list(knot_vectors_raw, tot_size)

        # Then reconstructor the knot vectors from the long vector by utilizing that we know how long each vector is
        knot_vectors = []
        start = 0
        for knot_vector_size in knot_vector_sizes:
            knot_vectors.append(knot_vectors_serialized[start:start+knot_vector_size])
            start += knot_vector_size

        return knot_vectors

    def get_control_points(self) -> ListList:
        """
        Get the matrix with the control points of the BSpline.
        :return Matrix (as a list of lists) with len(getCoefficients) rows and getNumOutputs columns
        """
        control_points_raw = splinter._call(splinter._get_handle().splinter_bspline_get_control_points, self._handle)

        num_rows = splinter._call(splinter._get_handle().splinter_bspline_get_num_control_points, self._handle)
        num_cols = self._num_variables + 1

        control_points_flattened = c_array_to_list(control_points_raw, num_rows * num_cols)

        control_points = []
        start = 0
        for row in range(num_rows):
            control_points.append(control_points_flattened[start:start+num_cols])
            start += num_cols

        return control_points

    def get_knot_averages(self) -> ListList:
        """
        Get the matrix with the knot averages of the BSpline.
        :return Matrix (as a list of lists) with len(getControlPoints) rows and getNumVariables columns
        """
        knot_averages_raw = splinter._call(splinter._get_handle().splinter_bspline_get_knot_averages, self._handle)

        num_rows = splinter._call(splinter._get_handle().splinter_bspline_get_num_control_points, self._handle)
        num_cols = self._num_variables

        knot_averages_flattened = c_array_to_list(knot_averages_raw, num_rows * num_cols)

        knot_averages = []
        start = 0
        for row in range(num_rows):
            knot_averages.append(knot_averages_flattened[start:start+num_cols])
            start += num_cols

        return knot_averages

    def get_basis_degrees(self) -> List[int]:
        """
        :return List with the basis degrees of the BSpline
        """
        num_vars = splinter._call(splinter._get_handle().splinter_bspline_get_dim_x, self._handle)
        basis_degrees = splinter._call(splinter._get_handle().splinter_bspline_get_basis_degrees, self._handle)

        return c_array_to_list(basis_degrees, num_vars)

    def insert_knots(self, val: float, dim: int, multiplicity: int=1):
        """
        Insert knot at 'val' to knot vector for variable 'dim'. The knot is inserted until a knot multiplicity of
        'multiplicity' is obtained.
        """
        splinter._call(splinter._get_handle().splinter_bspline_insert_knots, self._handle, val, dim, multiplicity)

    def decompose_to_bezier_form(self):
        """
        Insert knots until all knots have multiplicity degree + 1. This ensures that the polynomial pieces are not
        overlapping.
        """
        splinter._call(splinter._get_handle().splinter_bspline_decompose_to_bezier_form, self._handle)
