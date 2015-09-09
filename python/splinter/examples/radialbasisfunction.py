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

# Create a RadialBasisFunction of type Thin plate spline
rbf = splinter.RadialBasisFunction(d, splinter.RBFType.THIN_PLATE_SPLINE)

print("Jacobian at [3.2,3.2]: " + str(rbf.evalJacobian([3.2,3.2])))
print("Hessian at [3.2,3.2]: " + str(rbf.evalHessian([3.2,3.2])))

# Save the rbf to test.rbf
# The file ending doesn't matter
rbf.save("test.rbf")

# Create RadialBasisFunction from saved RadialBasisFunction
rbf2 = splinter.RadialBasisFunction("test.rbf")

print("Original RadialBasisFunction at [2.1,2.9]: " + str(rbf.eval([2.1,2.9])))
print("Loaded RadialBasisFunction at [2.1,2.9]: " + str(rbf2.eval([2.1,2.9])))

remove("test.rbf")