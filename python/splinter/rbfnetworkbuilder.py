# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from . import splinter
from .rbfnetwork import RBFNetwork
from .builderbase import BuilderBase


class RBFNetworkBuilder(BuilderBase):
    class RBFType:
        THIN_PLATE_SPLINE, MULTIQUADRIC, INVERSE_QUADRIC, INVERSE_MULTIQUADRIC, GAUSSIAN = range(1, 6)

        @staticmethod
        def isValid(value):
            return value in range(1, 6)

    def __init__(self, data):
        super(RBFNetworkBuilder, self).__init__(data)

        self._type = RBFNetworkBuilder.RBFType.THIN_PLATE_SPLINE
        self._normalized = False
        self._precondition = False

        self._handle = splinter._call(splinter._getHandle().splinter_rbfnetwork_builder_init, self._datatable._getHandle())

    def type(self, type):
        if not RBFNetworkBuilder.RBFType.isValid(type):
            raise ValueError("RBFNetworkBuilder:type: Invalid type. See RBFNetworkBuilder.RBFType for valid types.")

        self._type = type

        splinter._call(splinter._getHandle().splinter_rbfnetwork_builder_set_type, self._handle, self._type)
        return self

    def normalized(self, normalized):
        self._normalized = normalized

        splinter._call(splinter._getHandle().splinter_rbfnetwork_builder_set_normalized, self._handle, self._normalized)
        return self

    def precondition(self, precondition):
        self._precondition = precondition

        splinter._call(splinter._getHandle().splinter_rbfnetwork_builder_set_precondition, self._handle, self._precondition)
        return self

    # Returns a handle to the created internal RBFNetwork object
    def build(self):
        rbfnetwork_handle = splinter._call(splinter._getHandle().splinter_rbfnetwork_builder_build, self._handle)

        return RBFNetwork(rbfnetwork_handle)

    def __del__(self):
        if self._handle is not None:
            splinter._call(splinter._getHandle().splinter_rbfnetwork_builder_delete, self._handle)
        self._handle = None
