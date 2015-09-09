# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# Add the SPLINTER directory to the search path, so we can include it
from os import sys, path, remove
sys.path.append(path.dirname(path.dirname(path.dirname(path.abspath(__file__)))))


##### Start of the example #####
import splinter

splinter.load("/home/anders/SPLINTER/build/release/libsplinter-matlab-1-4.so")

def f(x):
	return x[0]*x[1]

# Create a DataTable and populate it with samples
d = splinter.DataTable()
for i in range(10):
	for j in range(10):
		d.addSample([i,j], f([i,j]))

# Create a BSpline of degree 1
b = splinter.BSpline(d, 1)

print("Jacobian at [3.2,3.2]: " + str(b.evalJacobian([[3.2,3.2],[1.0, 1.0]])))
print("Hessian at [3.2,3.2]: " + str(b.evalHessian([[3.2,3.2],[1.0, 1.0]])))

# Save the bspline to test.bspline
# The file ending doesn't matter
b.save("test.bspline")

# Create BSpline from saved BSpline
b2 = splinter.BSpline("test.bspline")

print("Original BSpline at [2.1,2.9]: " + str(b.eval([[2.1,2.9],[1.0,1.0],[2.0,2.0]])))
print("Loaded BSpline at [2.1,2.9]: " + str(b2.eval([[2.1,2.9]])))

remove("test.bspline")