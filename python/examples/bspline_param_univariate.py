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
import splinterpy

# Only for dev purposes
import os.path
if os.path.isdir("/home/bjarne/"):
    splinterpy.load("/home/bjarne/Code/splinter/build/debug/libsplinter-4-0.so")
elif os.path.isdir("/home/anders/"):
    splinterpy.load("/home/anders/SPLINTER/build/debug/libsplinter-4-0.so")

# B-spline built from parameters: coefficients, knot vectors and degrees
control_points = [0, 1, 0, 1, 0]
knot_vector = [0, 0, 1, 2, 3, 4, 4]
degree = 1
bs = splinterpy.BSpline.from_param(degree, knot_vector, control_points)

xd = np.arange(0, 4, .01)
yd = bs.eval(xd)

plt.plot(xd, yd, label='B-spline')
plt.legend(loc='upper right')
plt.show()
