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
        # If string we load the BSpline from the file
        if isString(handleOrFileName):
            fileName = getCString(handleOrFileName)
            self._handle = splinter._call(splinter._getHandle().bspline_load_init, fileName)

        # Else, the argument is the handle to the internal BSpline object
        else:
            self._handle = handleOrFileName

        self._numVariables = splinter._call(splinter._getHandle().function_get_num_variables, self._handle)
