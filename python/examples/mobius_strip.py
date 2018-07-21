# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Add the SPLINTER directory to the search path, so we can include it
import numpy as np
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinterpy
from math import sqrt, pi, cos, sin

# Only for dev purposes
import os.path
if os.path.isdir("/home/bjarne/"):
    splinterpy.load("/home/bjarne/Code/splinter/build/debug/libsplinter-4-0.so")
elif os.path.isdir("/home/anders/"):
    splinterpy.load("/home/anders/SPLINTER/build/debug/libsplinter-4-0.so")

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# Sample Möbius strip
u = np.linspace(0, 2*pi, 10)
v = np.linspace(-1, 1, 10)
n = len(u) * len(v)
X = np.ndarray((n, 2))
Y = np.ndarray((n, 3))
i = 0
for ui in u:
    for vi in v:
        X[i, :] = [ui, vi]
        Y[i, :] = [(1 + vi/2*cos(ui/2)) * cos(ui),
                   (1 + vi/2*cos(ui/2)) * sin(ui),
                   vi/2*sin(ui/2)]
        i = i + 1

# Build 3-D-valued bicubic B-spline with clamped knot vectors
bspline = splinterpy.bspline_interpolator(X.tolist(), Y.tolist(), degree=3)

# Sample B-spline and plot
u2 = np.linspace(0, 2*pi, 100)
v2 = np.linspace(-1, 1, 100)
U2, V2 = np.meshgrid(u2, v2)
x = np.ndarray((len(u2), len(v2)))
y = np.ndarray((len(u2), len(v2)))
z = np.ndarray((len(u2), len(v2)))

for i in range(len(u2)):
    for j in range(len(v2)):
        val = bspline.eval([U2[j, i], V2[j, i]])
        x[i, j] = val[0]
        y[i, j] = val[1]
        z[i, j] = val[2]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z, antialiased=True, )
plt.show()


