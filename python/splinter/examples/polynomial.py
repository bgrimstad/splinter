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
import numpy as np

# Must be done if splinter was unable to locate the shared library by itself
splinter.load("/home/anders/SPLINTER/build/release/libsplinter-2-0.so")


def f(x):
    return x[0]*x[1]

# Create a numpy array and populate it with samples
d = np.empty((10*10, 3))
idx = 0
for i in range(10):
    for j in range(10):
        d[idx] = [i, j, f([i, j])]
        idx += 1

# powers should be a matrix with the number of terms rows and number of variables columns
# If you want a Polynomial of the form f(x, y) = c1*x**2*y + c2*x*y**2, then powers should be [[2, 1], [1, 2]]
# This is because the first term is of degree 2 in x, degree 1 in y. The second term is of degree 1 in x, degree 2 in y.
powers = [1, 1]
poly = splinter.PolynomialBuilder(d).powers(powers).build()

print("Jacobian at [3.2,3.2]: " + str(poly.evalJacobian([3.2,3.2])))

# evalHessian is not implemented for Polynomial, expecting error:
try:
    print("Hessian at [3.2,3.2]: " + str(poly.evalHessian([3.2,3.2])))
except Exception as e:
    print(e)


# Save the Polynomial to test.poly
# The file ending doesn't matter
poly.save("test.poly")

# Create Polynomial from saved Polynomial
poly2 = splinter.Polynomial("test.poly")

print("Original Polynomial at [2.1,2.9]: " + str(poly.eval([2.1,2.9])))
print("Loaded Polynomial at [2.1,2.9]: " + str(poly2.eval([2.1,2.9])))

remove("test.poly")
