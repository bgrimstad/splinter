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
import os.path
if os.path.isdir("/home/bjarne/"):
    splinter.load("/home/bjarne/Code/C++/splinter/build/release/libsplinter-3-1.so")
elif os.path.isdir("/home/anders/"):
    splinter.load("/home/anders/SPLINTER/build/debug/libsplinter-3-1.so")


# Example showing a B-spline curve with a clamped knot vector
degree = [3]
knots = [[0, 0, 0, 0, 1, 1, 1, 1]]
control_points = [[0, 0], [1, 1], [2, 1], [3, 0]]
bspline_curve = splinter.BSpline.init_from_param(control_points, knots, degree)

# Evaluate B-spline curve for u in [0, 1]
u = np.linspace(0, 1, 1000)
y = [0] * len(u)
for i in range(len(u)):
    y[i] = bspline_curve.eval(u[i])

# Now y contains a list of 2-D coordinates
# The following call does not produce the same y!
y = bspline_curve.eval(u)

# Unzip to obtain two lists, one for each coordinate
y0, y1 = zip(*y)
p0, p1 = zip(*control_points)

# Plot results
plt.plot(p0, p1, '*', label='Control points')
plt.plot(y0, y1, '-', label='B-spline curve')
plt.legend(loc='upper left')
plt.show()
