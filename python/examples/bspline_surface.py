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


# Example showing a B-spline surface f : R^2 -> R^3 with clamped knot vectors
# f(u, v) = sum_i sum_j P_{ij} * B_i(u) * B_j(v),
# where P_{ij} is a control point in R^3
degrees = [3, 3]
knots = [[0, 0, 0, 0, 1, 1, 1, 1],
         [0, 0, 0, 0, 1, 1, 1, 1]]

control_points = []
for u in np.linspace(0, 1, len(knots[0]) - degrees[0] - 1):
    for v in np.linspace(0, 1, len(knots[1]) - degrees[1] - 1):
        control_points.append([u, v, sin(10*u) + 0.01*v])

bspline_surface = splinterpy.BSpline.from_param(degrees, knots, control_points)

# Evaluate
u = np.linspace(0, 1, 40)
v = np.linspace(0, 1, 40)
y = []

for u_i in u:
    for v_j in v:
        y.append(bspline_surface.eval([u_i, v_j]))

y0, y1, y2 = zip(*y)

# Plot B-spline surface
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_trisurf(y0, y1, y2, linewidth=0)

# Plot control structure
p0, p1, p2 = zip(*control_points)
p0_mesh, p1_mesh = np.meshgrid(p0, p1)
# ax.plot_wireframe(p0_mesh, p1_mesh, p2, color='black', alpha='0.2')
ax.scatter(p0, p1, p2)

plt.show()
