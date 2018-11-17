# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import copyreg
import tempfile

from .splinter_backend import splinter_backend_obj

from .bspline import BSpline
from .bsplineboosting import BSplineBoosting
from .bsplinebuilders import bspline_interpolator, bspline_smoother, bspline_unfitted


def load(lib_file_path: str):
    """
    Attempt to load the SPLINTER back-end from the file at lib_file_path.
    lib_file_path should be the shared library file (.so on Linux, etc.)
    :param lib_file_path:
    :return:
    """
    splinter_backend_obj.load(lib_file_path)


# Tell Pickle how to pickle and unpickle the BSpline class
def constructor(serialized_data) -> BSpline:
    with tempfile.NamedTemporaryFile() as temp:
        with open(temp.name, "wb") as f:
            f.write(serialized_data)

        return BSpline.from_json(temp.name)


def reducer(bspline: BSpline):
    with tempfile.NamedTemporaryFile() as temp:
        bspline.to_json(temp.name)
        with open(temp.name, "rb") as f:
            data = f.read()

    return constructor, (data, )


# Register reducer as the reducer for objects of type BSpline
copyreg.pickle(BSpline, reducer)

__all__ = [
    "bspline",
    "bsplinebuilders",
    "bsplineboosting",
]

try:
    splinter_backend_obj.load()
except Exception as e:
    print(e)

# Important to set this after the SPLINTER backend has been loaded
__version__ = splinter_backend_obj.version
