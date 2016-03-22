# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from os import sys, path, remove


# Add the SPLINTER directory to the search path, so we can include it
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

import splinter

# Load SPLINTER
splinter.load("/home/bjarne/Code/C++/splinter4/splinter/bin/Release/libsplinter-2-0.so")


# Example with two variables
def f(x):
    return x[0]*x[1]

x = np.arange(-5, 5, 0.25)
y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(x, y)
Z = np.sqrt(X**2)
Z = np.sqrt(X**2 + Y**2)

data = []
for i in range(len(x)):
    for j in range(len(y)):
        xij = X[i, j]
        yij = Y[i, j]
        zij = Z[i, j]
        data.append([xij, yij, zij])

# Cubic B-spline
bspline = splinter.BSplineBuilder(data, degree=[1, 3], smoothing=splinter.BSplineBuilder.Smoothing.NONE).build()

Zbs = Z

for i in range(len(x)):
    for j in range(len(y)):
        xij = X[i, j]
        yij = Y[i, j]
        Zbs[i, j] = bspline.eval([xij, yij])[0]

# Plot f
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=5)

# Plot b-spline
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Zbs, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=5)

# TODO: plot grid using plot_wireframe()

plt.show()

print("Jacobian at [3.2,3.2] and [1.0, 1.0]: " + str(bspline.evalJacobian([[3.2, 3.2], [1.0, 1.0]])))
print("Hessian at [3.2,3.2] and [1.0, 1.0]: " + str(bspline.evalHessian([[3.2, 3.2], [1.0, 1.0]])))

# Save the bspline to test.bspline
# The file ending doesn't matter
bspline.save("test.bspline")

# Create BSpline from saved BSpline
bspline = splinter.BSpline("test.bspline")

print("Original BSpline at [2.1,2.9],[1.0,1.0] and [2.0,2.0]: " + str(bspline.eval([[2.1, 2.9], [1.0, 1.0], [2.0, 2.0]])))
print("Loaded BSpline at [2.1,2.9],[1.0,1.0] and [2.0,2.0]: " + str(bspline.eval([[2.1, 2.9], [1.0, 1.0], [2.0, 2.0]])))

remove("test.bspline")