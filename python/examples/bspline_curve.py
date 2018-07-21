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


# Example showing a B-spline curve with a clamped knot vector
degree = 3
knots = [0, 0, 0, 0, 1, 1, 1, 1]
control_points = [[0, 0], [1, 1], [2, 1], [3, 0]]
bspline_curve = splinterpy.BSpline.from_param(degree, knots, control_points)

# Evaluate B-spline curve for u in [0, 1]
u = np.linspace(0, 1, 1000)
y = bspline_curve.eval(u)

# Now y contains a list of 2-D coordinates

# Unzip to obtain two lists, one for each coordinate
y0, y1 = zip(*y)
p0, p1 = zip(*control_points)

# Plot results
plt.plot(p0, p1, '*', label='Control points')
plt.plot(y0, y1, '-', label='B-spline curve')
plt.legend(loc='upper left')
plt.show()
