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
import pandas as pd

# Must be done if splinter was unable to locate the shared library by itself
splinter.load("/home/anders/SPLINTER/build/release/libsplinter-2-0.so")


def f(x0, x1):
    import math
    return x1 * (x0 - math.log(x0))

# Create a Pandas DataFrame and populate it with samples
d = []
for x0 in range(1, 10):
    for x1 in range(1, 10):
        d.append([x0, x1, f(x0, x1)])
df = pd.DataFrame(data=d, columns=["x0", "x1", "f(x0, x1)"])

# Create a RBFNetwork of type Thin plate spline
rbf = splinter.RBFNetworkBuilder(df.values).build()

print("Jacobian at [3.2, 3.2]: " + str(rbf.evalJacobian([3.2, 3.2])))

# evalHessian is not implemented for RBFNetwork, expecting error:
try:
    print("Hessian at [3.2, 3.2]: " + str(rbf.evalHessian([3.2, 3.2])))
except Exception as e:
    print(e)

# Save the rbf to test.rbf
# The file ending doesn't matter
rbf.save("test.rbf")

# Create RBFNetwork from saved RBFNetwork
rbf2 = splinter.RBFNetwork("test.rbf")

print("Original RBFNetwork at [2.1,2.9]: " + str(rbf.eval([2.1, 2.9])))
print("Loaded RBFNetwork at [2.1,2.9]: " + str(rbf2.eval([2.1, 2.9])))

remove("test.rbf")
