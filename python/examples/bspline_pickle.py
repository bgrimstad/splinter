# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import pickle

import splinterpy


filename = "test.pickle"

try:
    f = open(filename, "rb")
    bs = pickle.load(f)
    print(f"Loaded BSpline from {filename}")

except FileNotFoundError:
    f = open(filename, "wb")
    bs = splinterpy.bspline_interpolator([0, 1, 2], [2, 3, 4], 0)

    pickle.dump(bs, f)
    print(f"Saved BSpline to {filename}")

print(bs.eval(0))
