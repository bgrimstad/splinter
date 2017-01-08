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
import splinter
# Only for dev purposes
splinter.load("/home/bjarne/Code/C++/splinter/bin/Release/libsplinter-3-1.so")
# splinter.load("/home/anders/SPLINTER/build/debug/libsplinter-3-0.so")


def underlying_func(x: np.array):
    return 0.5*x + 0.1*np.apply_along_axis(np.square, 0, x)


def noisy_func(x: np.array):
    np.random.seed(123456)
    return underlying_func(x) + np.random.randn(x.shape[0])


print("Testing stochastic gradient boosting with B-spline learners")
fh = splinter.bsplineboosting.mse
a = np.array([3, 4, 5])
b = np.array([2, 3, 3])
print(fh(a, b))
print(splinter.bsplineboosting.ls(b))

# Sampling
x = np.arange(0, 10, 0.1)
xd = np.arange(-2, 12, 0.01)
y = noisy_func(x)
yd = underlying_func(xd)

# Just one P-spline
pspline = splinter.BSplineBuilder(x, y,
                                  smoothing=splinter.BSplineBuilder.Smoothing.PSPLINE,
                                  alpha=10.0,
                                  knot_spacing=splinter.BSplineBuilder.KnotSpacing.EXPERIMENTAL,
                                  num_basis_functions=20).build()

yd_pspline = pspline.eval(xd)

# Boosting
bb = splinter.BSplineBoosting(learning_rate=0.05, n_estimators=100, alpha=10.0)
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
