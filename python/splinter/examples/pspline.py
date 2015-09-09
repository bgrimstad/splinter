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

# Must be done if splinter was unable to locate the shared library by itself
splinter.load("/home/anders/SPLINTER/build/release/libsplinter-matlab-1-4.so")

def f(x):
	return x[0]*x[1]

# Create a DataTable and populate it with samples
d = splinter.DataTable()
for i in range(10):
	for j in range(10):
		d.addSample([i,j], f([i,j]))

# Create a PSpline with lambda (or smoothing parameter) 0.03
p = splinter.PSpline(d, 0.03)

print("Jacobian at [3.2,3.2]: " + str(p.evalJacobian([3.2,3.2])))
print("Hessian at [3.2,3.2]: " + str(p.evalHessian([3.2,3.2])))

# Save the pspline to test.pspline
# The file ending doesn't matter
p.save("test.pspline")

# Create PSpline from saved PSpline
p2 = splinter.PSpline("test.pspline")

print("Original PSpline at [2.1,2.9]: " + str(p.eval([2.1,2.9])))
print("Loaded PSpline at [2.1,2.9]: " + str(p2.eval([2.1,2.9])))

remove("test.pspline")