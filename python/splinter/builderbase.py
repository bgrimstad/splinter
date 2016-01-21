# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from ctypes import *
from .datatable import DataTable


class BuilderBase:
    def __init__(self, data):
        self._handle = None  # Handle for referencing the c side of this object

        self._datatable = DataTable(data)
        self._lambda = 0.03

    # These functions must be implemented by subclasses
    def setLambda(self, newLambda):
        pass

    def build(self):
        pass

    def __del__(self):
        pass
