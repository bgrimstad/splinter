# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import copyreg as _copyreg
import tempfile as _tempfile

from .splinter_backend import splinter_backend_obj as _backend

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
    _backend.load(lib_file_path)


# Tell Pickle how to pickle and unpickle the BSpline class
def _constructor(serialized_data) -> BSpline:
    with _tempfile.NamedTemporaryFile() as temp:
        with open(temp.name, "wb") as f:
            f.write(serialized_data)

        return BSpline.from_json(temp.name)


def _reducer(bspline: BSpline):
    with _tempfile.NamedTemporaryFile() as temp:
        bspline.to_json(temp.name)
        with open(temp.name, "rb") as f:
            data = f.read()

    return _constructor, (data,)


# Register reducer as the reducer for objects of type BSpline
_copyreg.pickle(BSpline, _reducer)

__all__ = [
    "bspline",
    "bsplinebuilders",
    "bsplineboosting",
]

try:
    _backend.load()
except Exception as e:
    print(e)

# Important to set this after the SPLINTER backend has been loaded
__version__ = _backend.version
