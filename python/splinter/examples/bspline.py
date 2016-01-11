# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# Add the SPLINTER directory to the search path, so we can include it
import numpy as np
import matplotlib.pyplot as plt
from os import sys, path, remove
sys.path.append(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))


##### Start of the example #####
import splinter

splinter.load("/home/bjarne/Code/C++/splinter4/splinter/bin/Release/libsplinter-2-0.so")


# Example with one variables
def f1(x):
    return -1. + 2*x + 0.1*(x**2) + 10*np.random.rand(1)[0]

# Create a DataTable and populate it with samples
x = [0.]*11
y = [0.]*11
d1 = splinter.DataTable()
for i in range(11):
    x[i] = i
    y[i] = f1(i)
    d1.addSample(x[i], y[i])

# Cubic B-spline that interpolates the data
b1 = splinter.BSplineBuilder(d1)\
    .smoothing(splinter.BSplineBuilder.Smoothing.NONE)\
    .build()

# .degree([splinter.BSplineBuilder.Degree.LINEAR, splinter.BSplineBuilder.Degree.CUBIC])\

# Cubic B-spline with regularization
b2 = splinter.BSplineBuilder(d1)\
    .smoothing(splinter.BSplineBuilder.Smoothing.REGULARIZATION)\
    .setLambda(0.1)\
    .build()

# Cubic P-spline
b3 = splinter.BSplineBuilder(d1)\
    .smoothing(splinter.BSplineBuilder.Smoothing.PSPLINE)\
    .setLambda(0.1)\
    .build()

n = 1000
xd = [0.]*n
yd1 = [0.]*n
yd2 = [0.]*n
yd3 = [0.]*n
for i in range(n):
    val = i/100.
    xd[i] = val
    yd1[i] = b1.eval([val])
    yd2[i] = b2.eval([val])
    yd3[i] = b3.eval([val])


plt.plot(x, y, '*')
plt.plot(xd, yd1)
plt.plot(xd, yd2, '--')
plt.plot(xd, yd3, '-.')
plt.show()


# Example with two variables
def f2(x):
    return x[0]*x[1]

# Create a DataTable and populate it with samples
d2 = splinter.DataTable()
for i in range(11):
    for j in range(11):
        d2.addSample([i,j], f2([i,j]))

# Cubic P-spline
b2d = splinter.BSplineBuilder(d2)\
    .smoothing(splinter.BSplineBuilder.Smoothing.NONE)\
    .build()

print("Jacobian at [3.2,3.2] and [1.0, 1.0]: " + str(b2d.evalJacobian([[3.2,3.2], [1.0, 1.0]])))
print("Hessian at [3.2,3.2] and [1.0, 1.0]: " + str(b2d.evalHessian([[3.2,3.2], [1.0, 1.0]])))

# Save the bspline to test.bspline
# The file ending doesn't matter
b2d.save("test.bspline")

# Create BSpline from saved BSpline
b2d = splinter.BSpline("test.bspline")

print("Original BSpline at [2.1,2.9],[1.0,1.0] and [2.0,2.0]: " + str(b2d.eval([[2.1,2.9], [1.0,1.0], [2.0,2.0]])))
print("Loaded BSpline at [2.1,2.9],[1.0,1.0] and [2.0,2.0]: " + str(b2d.eval([[2.1,2.9], [1.0,1.0], [2.0,2.0]])))

remove("test.bspline")