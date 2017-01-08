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
import splinter

# Only for dev purposes
splinter.load("/home/bjarne/Code/C++/splinter/bin/Release/libsplinter-3-1.so")
# splinter.load("/home/anders/SPLINTER/build/debug/libsplinter-3-0.so")

# B-spline built from parameters: coefficients, knot vectors and degrees
single_knot_vector = [k for k in np.linspace(-10, 10, 21)]
knot_vectors = [single_knot_vector] * 3

degrees = [1, 2, 3]

num_control_points = 1
for i, kv in enumerate(knot_vectors):
    num_control_points *= len(kv) - degrees[i] - 1

control_points = np.linspace(0, 100, num_control_points)

bs = splinter.BSpline.init_from_param(control_points, knot_vectors, degrees)

xd = [0, 0, 0]
yd = bs.eval(xd)
print(yd)

knot_averages = bs.get_knot_averages()
control_points = bs.get_control_points()
print(knot_averages)
