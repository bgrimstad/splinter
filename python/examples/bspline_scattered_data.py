# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Add the SPLINTER directory to the search path, so we can include it
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import sin
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinterpy

# Only for dev purposes
import os.path
if os.path.isdir("/home/bjarne/"):
    splinterpy.load("/home/bjarne/Code/splinter/build/debug/libsplinter-4-0.so")
elif os.path.isdir("/home/anders/"):
    splinterpy.load("/home/anders/SPLINTER/build/debug/libsplinter-4-0.so")


# Example showing how to fit a B-spline f : D -> R to scattered data on a domain D in R^2

# First let us define a zero-valued B-spline on D = [0, 5] x [0, 5]
degrees = [3, 3]
knots = [[0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5],
         [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]]
dim_y = 1
bspline = splinterpy.BSpline.from_param(degrees, knots, dim_y)

# Draw some samples randomly on D
num_samples = 10
np.random.seed(3)
x = np.random.rand(num_samples, 2) * 5  # Draw samples (x0, x1) in D
x = x.tolist()
y = np.random.rand(num_samples)  # Draw samples in [0, 1)
weights = [1.0] * num_samples
fitted_bspline = bspline.fit(x, y)

# Evaluate fitted B-spline
x0_e = []
x1_e = []
y_e = []

for u in np.linspace(0, 5, 40):
    for v in np.linspace(0, 5, 40):
        x0_e.append(u)
        x1_e.append(v)
        y_e.append(fitted_bspline.eval([u, v]))

# Plot B-spline surface
x0 = [x0 for x0, x1 in x]
x1 = [x1 for x0, x1 in x]
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_trisurf(x0_e, x1_e, y_e, linewidth=0, alpha=0.5)
ax.scatter(x0, x1, y, marker='x', s=40, c='black', alpha=1.0)
plt.show()
