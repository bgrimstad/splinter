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
import splinter_py

# Only for dev purposes
import os.path
if os.path.isdir("/home/bjarne/"):
    splinter_py.load("/home/bjarne/Code/C++/splinter/build/release/libsplinter-4-0.so")
elif os.path.isdir("/home/anders/"):
    splinter_py.load("/home/anders/SPLINTER/build/debug/libsplinter-4-0.so")


# Example with one variable
def f1(x):
    return -1. + 2*x + 0.1*(x**2) + 10*np.random.rand(1)[0]

x = np.arange(0, 11, 1)
y = np.array([f1(x_i) for x_i in x])

# Piecewise constant B-spline that interpolates the data
b0 = splinter_py.BSplineBuilder(1, 1, degree=0).fit(x, y)
print(b0.get_control_points())
print(b0.get_knot_vectors())

# Linear B-spline that interpolates the data
b1 = splinter_py.BSplineBuilder(1, 1, degree=1).fit(x, y)

# Quadratic B-spline that interpolates the data
b2 = splinter_py.BSplineBuilder(1, 1, degree=2).fit(x, y)

# Cubic B-spline that interpolates the data
b3 = splinter_py.BSplineBuilder(1, 1, degree=3).fit(x, y)

xd = np.arange(0, 10, .01)
yd0 = b0.eval(xd)
yd1 = b1.eval(xd)
yd2 = b2.eval(xd)
yd3 = b3.eval(xd)

plt.plot(x, y, '*', label='Data points')
plt.plot(xd, yd0, label='Piecewise constant B-spline')
plt.plot(xd, yd1, label='Linear B-spline')
plt.plot(xd, yd2, '--', label='Quadratic B-spline')
plt.plot(xd, yd3, '-.', label='Cubic B-spline')
plt.legend(loc='upper left')
plt.show()
