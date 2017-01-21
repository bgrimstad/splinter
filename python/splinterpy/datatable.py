# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from .splinter_backend import splinter_backend_obj
from .utilities import *


class DataTable:
    """
    This class is only intended for internal use by the SPLINTER Python interface
    """
    def __init__(self, x_or_data, y=None):
        self.__handle = None  # Handle to instance in the library
        self.__dim_x = None
        self.__dim_y = None
        self.__num_samples = 0  # Number of samples not yet transferred to back end
        self.__xs = []
        self.__ys = []

        if is_string(x_or_data):
            self.__handle = splinter_backend_obj.call(splinter_backend_obj.handle.splinter_datatable_load_init, get_c_string(x_or_data))
            self.__dim_x = splinter_backend_obj.call(splinter_backend_obj.handle.splinter_datatable_get_dim_x, self.__handle)
            self.__dim_y = splinter_backend_obj.call(splinter_backend_obj.handle.splinter_datatable_get_dim_y, self.__handle)
        else:
            self.__handle = splinter_backend_obj.call(splinter_backend_obj.handle.splinter_datatable_init)

            if y is None:
                raise Exception("No y-values supplied.")

            if len(x_or_data) != len(y):
                raise Exception("x and y must be of the same length!")

            # If x_or_data is not a string, we expect it to be lists of x values which has corresponding y values in 'y'.
            for x, y in zip(x_or_data, y):
                self.add_sample(x, y)

    def add_sample(self, x, y):
        if not isinstance(x, list):
            x = [x]

        if not isinstance(y, list):
            y = [y]

        # New DataTable, this sample is the first. Let it determine the dimensionality of the domain and codomain
        if self.__dim_x is None:
            self.__dim_x = len(x)
            self.__dim_y = len(y)

        if self.__dim_x != len(x):
            raise Exception("Dimension of the domain of the new sample disagrees with the dimension of the domain of previous samples!\n"
                            "Previous: {}, new: {}".format(self.__dim_x, len(x)))

        if self.__dim_y != len(y):
            raise Exception("Dimension of the codomain of the new sample disagrees with the dimension of the codomain of previous samples!\n"
                            "Previous: {}, new: {}".format(self.__dim_y, len(y)))

        x_float = [float(_x) for _x in x]
        y_float = [float(_y) for _y in y]

        self.__xs += x_float
        self.__ys += y_float
        self.__num_samples += 1

    def get_dim_x(self):
        return self.__dim_x

    def get_dim_y(self):
        return self.__dim_y

    def get_num_samples(self):
        self.__transfer()
        return splinter_backend_obj.call(splinter_backend_obj.handle.splinter_datatable_get_num_samples, self.__handle)

    # Getter for the datatable for use by BSpline
    # Will make sure all samples are transferred to the back end before returning the handle to a BSpline
    def _get_handle(self):
        self.__transfer()
        return self.__handle

    # Methods below are internal use by the class only
    # Transfer samples to the library
    def __transfer(self):
        if self.__num_samples > 0:
            func_handle = splinter_backend_obj.handle.splinter_datatable_add_samples_row_major
            splinter_backend_obj.call(func_handle, self.__handle,
                                      list_to_c_array_of_doubles(self.__xs), self.__dim_x,
                                      list_to_c_array_of_doubles(self.__ys), self.__dim_y,
                                      self.__num_samples)

            self.__num_samples = 0
