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


# B-spline built from parameters: coefficients, knot vectors and degrees
control_points = [0, 1, 0, 1, 0]
knot_vector = [0, 0, 1, 2, 3, 4, 4]
degree = 1

# control_points = [0, 1, 0]
# knot_vector = [0, 1, 2, 3, 4, 5]
# degree = 2

control_points = [0, 1, 0, 0, 1, 0]
knot_vector = [0, 0, 1, 2, 2, 3, 4, 4]
degree = 1

control_points = [0, 1, 0, 0, 1, 0]
knot_vector = [0, 0, 0, 1, 2, 2, 3, 4, 4, 4]
degree = 3

bs = splinterpy.BSpline.from_param(degree, knot_vector, control_points)


# Compute B-spline derivative (see pp. 93-94 in the NURBS book)
def compute_control_points(c, p, t):
    """
    Return unmasked control points of derivative
    The number of new control points is len(c) + 1
    :param c: control points
    :param p: degree
    :param t: knot vector
    :return: new control points
    """
    c_new = []

    # First control point
    dt = t[p] - t[0]
    if dt == 0:
        c_new.append(0)
    else:
        c_new.append(p*c[0]/dt)

    # Interior control points (Qi)
    for i in range(len(c) - 1):
        dt = t[i+p+1] - t[i+1]
        if dt == 0:
            c_new.append(0)
        else:
            c_new.append(p*(c[i+1] - c[i])/(t[i+p+1] - t[i+1]))

    # End control point
    n = len(c)
    dt = t[n+p] - t[n]
    if dt == 0:
        c_new.append(0)
    else:
        c_new.append(-p*c[-1]/dt)

    return c_new


# For all knots of multiplicity p+1, remove the first knot, bringing the multiplicity to p
def remove_multiples(knots, p):
    if len(knots) == 0:
        return []
    knots_new = [knots[0]]
    t_prev = knots[0]
    m_count = 1
    for t in knots[1:]:
        if t == t_prev:
            m_count += 1
        else:
            m_count = 1
        t_prev = t

        if m_count >= p + 1:
            continue
        else:
            knots_new.append(t)
    return knots_new


def control_point_mask(knots, p):
    """
    Assume that knot vector is sorted (regular)
    NOTE: a zero may not occur among the last p elements of the mask
    :param knots:
    :param p:
    :return: mask
    """
    cp_mask = []  # Mask out control points related to p+1 multiple knots
    unique_knots, knot_counts = np.unique(knots, return_counts=True)

    for unique, count in zip(unique_knots, knot_counts):
        if count > p + 1:
            raise ValueError("Knot vector not regular")
        elif count == p + 1:
            cp_mask += [0] + [1]*p
        else:
            cp_mask += [1]*count

    return cp_mask


def mask_control_points(c, mask):
    masked_control_points = []
    for i, ci in enumerate(c):
        if not mask[i] == 0:
            masked_control_points.append(ci)
    return masked_control_points


knot_vector2 = remove_multiples(knot_vector, degree)
print("New knot vector", knot_vector2)

knots_removed = len(knot_vector) - len(knot_vector2)
print("Knots removed:", knots_removed)

control_points2 = compute_control_points(control_points, degree, knot_vector)
print("New control points:", control_points2)

mask = control_point_mask(knot_vector, degree)
assert(len(mask) == len(knot_vector))
print("Mask:", mask)

masked_control_points = mask_control_points(control_points2, mask)
print("Masked control points:", masked_control_points)

# Build derivative
bs2 = splinterpy.BSpline.from_param(degree - 1, knot_vector2, masked_control_points)

# Evaluate and plot results
xd = np.arange(0, 4, .01)
yd = bs.eval(xd)

xd2 = np.arange(0, 4, .01)
yd2 = bs2.eval(xd2)

xd3 = np.arange(0, 4, .01)
yd3 = bs.eval_jacobian(xd3)

plt.plot(xd, yd, label='B-spline')
plt.legend(loc='upper right')

plt.plot(xd2, yd2, label='B-spline derivative')
plt.legend(loc='upper right')

plt.plot(xd3, yd3, '-*', label='B-spline pointwise derivative')
plt.legend(loc='upper right')

plt.show()
