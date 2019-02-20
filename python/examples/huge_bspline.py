# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import splinterpy as spy


np.random.seed(1234)

# Input dimension
dim = 4

# Create grid
x_sloc = 20  # Number of sample locations in each variable
x_start = 1
x_stop = 50
x = np.linspace(x_start, x_stop, x_sloc).tolist()

grid = [x] * dim
sloc = list(itertools.product(*grid))  # Sample locations (as list of tuples)
X = [list(s) for s in sloc]  # Convert to list since SPLINTER will complain about tuples

# Number of sample locations should be x_sloc**dim
print("Number of sample locations:", len(X))

# Draw random samples from normal distribution
y = np.random.normal(0, 1, len(X)).tolist()

# Build B-spline
bspline = spy.bspline_interpolator(X, y, degree=3)

# Evaluate
print(bspline.eval([1, 2, 3, 4]))
