# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Add the SPLINTER directory to the search path, so we can include it
import numpy as np
import matplotlib.pyplot as plt
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinter

# Only for dev purposes
splinter.load("/home/bjarne/Code/C++/splinter/splinter/bin/Release/libsplinter-3-1.so")
# splinter.load("/home/anders/SPLINTER/build/debug/libsplinter-3-0.so")

# B-spline built from parameters: coefficients, knot vectors and degrees
coefficients = [0, 1, 0, 1, 0]
knot_vectors = [[0, 0, 1, 2, 3, 4, 4]]
degrees = [1]
bs = splinter.BSpline.init_from_param(coefficients, knot_vectors, degrees)

xd = np.arange(0, 4, .01)
yd = bs.eval(xd)

plt.plot(xd, yd, label='B-spline')
plt.legend(loc='upper right')
plt.show()
