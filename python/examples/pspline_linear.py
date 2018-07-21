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


# Example with one variable
x = np.linspace(0, 10, 11)
y = 1*x + 1

# Piecewise constant B-spline that interpolates the data
PSPLINE = splinterpy.BSpline.Smoothing.PSPLINE
bs = splinterpy.bspline_smoother(x, y, degree=3, smoothing=PSPLINE, alpha=1.0)

xd = np.arange(0, 10, .01)
yd = bs.eval(xd)

plt.plot(x, y, '*', label='Data points')
plt.plot(xd, yd, label='P-spline')
plt.legend(loc='upper left')
plt.show()
