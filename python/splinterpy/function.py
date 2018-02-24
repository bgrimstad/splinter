# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from .splinter_backend import splinter_backend_obj
from .utilities import *
import numpy as np


class Function(object):
    def __init__(self):
        self._handle = None
        self._dim_x = None
        self._dim_y = None

    def eval(self, x):
        x = self._transform_input(x)

        num_points = len(x) // self._dim_x

        if len(x) % self._dim_x != 0:
            raise Exception("Attempting to evaluate function with input that is not a multiple of the dimension of the"
                            "function! {} % {} != 0".format(len(x), self._dim_x))

        f_handle = splinter_backend_obj.handle.splinter_bspline_eval_row_major
        res = splinter_backend_obj.call(f_handle, self._handle, list_to_c_array_of_doubles(x), len(x))

        results = []
        # Convert from ctypes array to Python list of lists
        # results is a list of the results in all evaluated points
        if self._dim_y == 1:
            results = c_array_to_list(res, num_points)
        else:
            for i in range(0, num_points*self._dim_y, self._dim_y):
                this_result = res[i:i+self._dim_y]
                results.append(c_array_to_list(this_result, len(this_result)))

        # We don't want to return a list when only evaluating one point
        if num_points == 1:
            return results[0]
        return results

    def eval_jacobian(self, x):
        x = self._transform_input(x)

        num_points = len(x) // self._dim_x

        f_handle = splinter_backend_obj.handle.splinter_bspline_eval_jacobian_row_major
        jac = splinter_backend_obj.call(f_handle, self._handle, list_to_c_array_of_doubles(x), len(x))

        # TODO: add support for self._dim_y > 1
        # Convert from ctypes array to Python list of lists
        # jacobians is a list of the jacobians in all evaluated points
        jacobians = []
        for i in range(num_points):
            jacobians.append([])
            for j in range(self._dim_x):
                jacobians[i].append(jac[i * self._dim_x + j])
        return jacobians

    def eval_hessian(self, x):
        x = self._transform_input(x)

        num_points = len(x) // self._dim_x

        f_handle = splinter_backend_obj.handle.splinter_bspline_eval_hessian_row_major
        hes = splinter_backend_obj.call(f_handle, self._handle, list_to_c_array_of_doubles(x), len(x))

        # Convert from ctypes array to Python list of list of lists
        # hessians is a list of the hessians in all points
        hessians = []
        for i in range(num_points):
            hessians.append([])
            for j in range(self._dim_x):
                hessians[i].append([])
                for k in range(self._dim_x):
                    hessians[i][j].append(hes[i * self._dim_x * self._dim_x + j * self._dim_x + k])
        return hessians

    def get_dim_x(self):
        return self._dim_x

    def get_dim_y(self):
        return self._dim_y

    @staticmethod
    def _transform_input(x):
        if isinstance(x, np.ndarray):
            x = x.tolist()

        if not isinstance(x, list):
            x = [x]

        # See if x is on the form [[x0,x1],[x2,x3]], if not we assume it to be on the form [x0, x1, x2, x3]
        if isinstance(x[0], list):
            x = flatten_list(x)

        return x

    def __del__(self):
        if self._handle is not None:
            splinter_backend_obj.call(splinter_backend_obj.handle.splinter_bspline_delete, self._handle)
        self._handle = None
