# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# This example uses the numpy and pillow libraries
# Install them using
# pip install numpy pillow
#
# This example loads an image, removes every nth row and column from it and fits three B-splines to the remaining data.
# The B-splines are then used to reconstruct the image in the original resolution, filling in the rows and columns we
# removed.
#
# A smoothing B-spline (penalized spline) could have been used as well, but without careful tuning of the
# alpha-parameter the results tend to get blurry due to the smoothing.

import splinterpy as spy

import numpy as np
from PIL import Image

# Start settings
original_img = Image.open("lenna.png")
# Keep every nth row and column from the image
decimation_factor = 4
# Use a cubic BSpline interpolator for reconstruction. degree can be any of 0, 1, 2, 3 and 4
degree = 3
# End settings


original = np.asarray(original_img)

original_dim_y, original_dim_x, _ = original.shape

# Create a fully black image of the same size as the original
reduced = np.zeros(original.shape, original.dtype)

# Fill reduced with the rows and columns from the original that we want to keep. The rest is kept black.
reduced[::decimation_factor, ::decimation_factor, :] = original[::decimation_factor, ::decimation_factor, :]

# Create a tensor containing only the data from the rows and columns that we are keeping for fitting the BSpline
ys = original[::decimation_factor, ::decimation_factor, :]


# Sample the original image at an equidistant grid with `decimation_factor` distance between sample points
xs_sample = []
for i in range(0, original_dim_y, decimation_factor):
    for j in range(0, original_dim_x, decimation_factor):
        xs_sample.append([i, j])

# Three color channels, we interpolate each separately from the others
color_channels = []
for i in range(3):
    ys_this_channel = ys[:, :, i].reshape((len(xs_sample),))

    # Try to change the call to bspline_interpolator to bspline_smoother
    interpolator = spy.bspline_interpolator(xs_sample, ys_this_channel, degree=degree)

    color_channels.append(interpolator)

# Create a list of all the points we want to sample the BSpline at
xs_eval = []
for i in range(original_dim_y):
    for j in range(original_dim_x):
        xs_eval.append([i, j])

# Evaluate the three BSplines in each of the points in xs_eval
restored_channels = []
for color_channel in color_channels:
    channel_values = color_channel.eval(xs_eval)
    channel_values = np.array(channel_values).reshape((original_dim_y, original_dim_x, 1))

    # The result of the evaluation may have values outside the valid range, which is [0, 255]
    channel_values = np.clip(channel_values, 0, 255)

    restored_channels.append(channel_values)

# Concatenate the channels to get a tensor of original_dim_x, original_dim_y, 3 shape
restored = np.concatenate(restored_channels, axis=2).astype(np.uint8)

# Show the original image, decimated image and restored image together for easy comparison
buf = np.zeros((original_dim_y, 3*original_dim_x, 3), dtype=np.uint8)
buf[:, :original_dim_x, :] = original
buf[:, original_dim_x:2*original_dim_x, :] = reduced
buf[:, 2*original_dim_x:3*original_dim_x, :] = restored

Image.fromarray(buf).show()
