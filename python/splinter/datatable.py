# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from .utilities import *


class DataTable:
    def __init__(self, filename_or_data=None):
        self.__handle = None  # Handle to instance in the library
        self.__x_dim = None
        self.__num_samples = 0  # Number of samples not yet transferred to back end
        self.__samples = []

        if is_string(filename_or_data):
            self.__handle = splinter._call(splinter._get_handle().splinter_datatable_load_init, get_c_string(filename_or_data))
            self.__x_dim = splinter._call(splinter._get_handle().splinter_datatable_get_num_variables, self.__handle)
        else:
            self.__handle = splinter._call(splinter._get_handle().splinter_datatable_init)

            # If fileNameOrData is not a string and not none, we expect it to be an iterable with rows of data,
            # where the last column in each row is y and the first n-1 columns is x.
            if filename_or_data is not None:
                for data_row in filename_or_data:
                    self.add_sample(list(data_row[:-1]), data_row[-1])

    # "Public" methods (for use by end user of the library)
    def add_sample(self, x, y):
        if not isinstance(x, list):
            x = [x]

        if self.__x_dim is None:
            self.__x_dim = len(x)

        if self.__x_dim != len(x):
            raise Exception("Dimension of the new sample disagrees with the dimension of previous samples!\nPrevious: " + str(self.__x_dim) + ", new: " + str(len(x)))

        self.__samples += list(x)
        self.__samples += [y]
        self.__num_samples += 1

    def get_num_variables(self):
        return self.__x_dim

    def get_num_samples(self):
        self.__transfer()
        return splinter._call(splinter._get_handle().splinter_datatable_get_num_samples, self.__handle)

    # Methods below are internal use only

    # Transfer samples to the library
    def __transfer(self):
        #print("Transferring " + str(self.__numSamples) + " samples to backend:")
        #for i in range(self.__numSamples):
        #	print(str(self.__samples[i*(self.__xDim+1)]) + "," + str(self.__samples[i*(self.__xDim+1)+1]) + " = " + str(self.__samples[i*(self.__xDim+1)+2]))

        if self.__num_samples > 0:
            splinter._call(splinter._get_handle().splinter_datatable_add_samples_row_major, self.__handle, (c_double * len(self.__samples))(*self.__samples), self.__num_samples, self.__x_dim)

        self.__samples = []
        self.__num_samples = 0

    # Getter for the datatable for use by BSpline
    # Will make sure all samples are transferred to the back end before returning the handle to a BSpline
    def _get_handle(self):
        self.__transfer()
        return self.__handle
