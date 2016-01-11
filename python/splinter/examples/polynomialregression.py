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
splinter.load("/home/anders/SPLINTER/build/debug/libsplinter-2-0.so")


def f(x):
    return x[0]*x[1]

# Create a DataTable and populate it with samples
d = splinter.DataTable()
for i in range(10):
    for j in range(10):
        d.addSample([i,j], f([i,j]))

# Create a PolynomialRegression of degree 1 in all dimensions
# You can also specify the degree of all dimensions by providing a list of numbers
# (preferably ints, as they will get casted to that anyway)
# The list has to be of length equal to the number of dimensions (the number returned by DataTable.getNumVariables())
# Ex: poly = splinter.PolynomialRegression(d, [1,2])
# Will make a PolynomialRegression of degree 1 in dimension 1, and degree 2 in dimension 2.
poly = splinter.PolynomialRegression(d, 1)

print("Jacobian at [3.2,3.2]: " + str(poly.evalJacobian([3.2,3.2])))
print("Hessian at [3.2,3.2]: " + str(poly.evalHessian([3.2,3.2])))

# Save the PolynomialRegression to test.poly
# The file ending doesn't matter
poly.save("test.poly")

# Create PolynomialRegression from saved PolynomialRegression
poly2 = splinter.PolynomialRegression("test.poly")

print("Original PolynomialRegression at [2.1,2.9]: " + str(poly.eval([2.1,2.9])))
print("Loaded PolynomialRegression at [2.1,2.9]: " + str(poly2.eval([2.1,2.9])))

remove("test.poly")