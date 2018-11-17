# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from .datatable import DataTable
from .bspline import BSpline
from .bsplinebuilders import bspline_interpolator, bspline_smoother, bspline_unfitted
from .bsplineboosting import BSplineBoosting

from .splinter_backend import splinter_backend_obj

try:
    splinter_backend_obj.load()
except Exception as e:
    print(e)

__version__ = splinter_backend_obj.__version__


__all__ = [
    "splinter_backend",
    "bspline",
    "bsplinebuilders",
    "bsplineboosting"
]


def load(lib_file_path: str):
    """
    Attempt to load the SPLINTER back-end from the file at lib_file_path.
    lib_file_path should be the shared library file (.so on Linux, etc.)
    :param lib_file_path:
    :return:
    """
    splinter_backend_obj.load(lib_file_path)
