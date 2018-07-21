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

# Drawing of Picasso's minimalistic dog using Bezier curves
# Inspired by: https://jeremykun.com/2013/05/11/bezier-curves-and-picasso/
control_points = [
    [180, 280], [183, 268], [186, 256], [189, 244],  # Front leg
    [191, 244], [290, 244], [300, 230], [339, 245],  # Tummy
    [340, 246], [350, 290], [360, 300], [355, 210],  # Back leg
    [353, 210], [370, 207], [380, 196], [375, 193],  # Tail
    [375, 193], [310, 220], [190, 220], [164, 205],  # Back
    [164, 205], [135, 194], [135, 265], [153, 275],  # Ear start
    [153, 275], [168, 275], [170, 180], [150, 190],  # Ear end + head
    [149, 190], [122, 214], [142, 204], [85, 240],   # Nose bridge
    [86, 240], [100, 247], [125, 233], [140, 238]    # Mouth
]

# Nine Bezier segments (9 knot intervals)
knot_vector = [0, 0, 0, 0,
               1, 1, 1, 1,
               2, 2, 2, 2,
               3, 3, 3, 3,
               4, 4, 4, 4,
               5, 5, 5, 5,
               6, 6, 6, 6,
               7, 7, 7, 7,
               8, 8, 8, 8,
               9, 9, 9, 9]
degree = 3
bs = splinterpy.BSpline.from_param(degree, knot_vector, control_points)

u = np.arange(knot_vector[0], knot_vector[-1], .01)
p = bs.eval(u)
x = [pi[0] for pi in p]
y = [-pi[1] for pi in p]  # y coordinates are given top-down

plt.figure(figsize=(10, 3))
plt.plot(x, y, color='black', linewidth=2, label="Picasso's dog")
plt.legend(loc='upper center')
# plt.savefig('picasso.pdf', format='pdf', dpi='300')
plt.show()

