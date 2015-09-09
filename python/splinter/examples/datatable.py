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

# This must be done if splinter could not locate the shared library
splinter.load("/home/anders/SPLINTER/build/release/libsplinter-matlab-1-4.so")

# Create a new, empty DataTable
d = splinter.DataTable()

# Populate it with samples
for i in range(10):
	for j in range(10):
		d.addSample([i,j], i*j)

# Save the samples for use later
d.save("test.datatable")

# Create DataTable from saved
d2 = splinter.DataTable("test.datatable")

print("Original:")
print("Dimension of the samples: " + str(d.getNumVariables()))
print("Number of samples in the library: " + str(d.getNumSamples()))

print("Loaded:")
print("Dimension of the samples: " + str(d2.getNumVariables()))
print("Number of samples in the library: " + str(d2.getNumSamples()))

remove("test.datatable")