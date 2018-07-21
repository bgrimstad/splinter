# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

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


def underlying_func(x: np.array):
    return 1.0 + np.apply_along_axis(np.sin, 0, x)


def noisy_func(x: np.array):
    np.random.seed(123456)
    return underlying_func(x) + np.random.randn(x.shape[0])


print("Testing stochastic gradient boosting with B-spline learners")
fh = splinterpy.bsplineboosting.mse
a = np.array([3, 4, 5])
b = np.array([2, 3, 3])
print(fh(a, b))
print(splinterpy.bsplineboosting.ls(b))

# Sampling
x = np.arange(0, 20, 0.1)
y = noisy_func(x)
xd = np.arange(-2, 22, 0.01)
yd = underlying_func(xd)

# Just one P-spline
EQUI_KNOTS = splinterpy.BSpline.KnotSpacing.EXPERIMENTAL
PSPLINE = splinterpy.BSpline.Smoothing.PSPLINE
pspline = splinterpy.bspline_unfitted(x, y, degrees=[3], knot_spacing=EQUI_KNOTS, num_basis_functions=[20])\
    .fit(x, y, smoothing=PSPLINE, alpha=10.0)

yd_pspline = pspline.eval(xd)

# Boosting
bb = splinterpy.BSplineBoosting(learning_rate=0.05, n_estimators=100, alpha=10.0)
bb.fit(x, y)

# Prediction
yd_boost = bb.eval(xd)

# Plotting
plt.plot(x, y, '*', label='Data points')
plt.plot(xd, yd, label='Unknown function')
plt.plot(xd, yd_pspline, label='P-spline')
plt.plot(xd, yd_boost, label='Boosted spline')
plt.legend(loc='upper left')
plt.show()
