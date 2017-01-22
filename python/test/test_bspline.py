import pytest
from os import sys, path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinterpy

# Only for dev purposes
import os.path
if os.path.isdir("/home/bjarne/"):
    splinterpy.load("/home/bjarne/Code/C++/splinter/build/release/libsplinter-4-0.so")
elif os.path.isdir("/home/anders/"):
    splinterpy.load("/home/anders/SPLINTER/build/debug/libsplinter-4-0.so")


test_param = [
    ([0, 1, 0, 1, 0], [0, 0, 1, 2, 3, 4, 4], 1, None),
    ([0, 1, 0, 1, 0], [0, 0, 1, 2, 3, 4, 4], [1], None),
    ([0, 1, 0, 1, 0], [[0, 0, 1, 2, 3, 4, 4]], 1, None),
    ([0, 1, 0, 1, 0], [[0, 0, 1, 2, 3, 4, 4]], 1, None),
    ([[0, 1, 0, 1, 0]], [0, 0, 1, 2, 3, 4, 4], 1, Exception),
    ([0, 1, 0, 1, 0, 1], [0, 0, 0, 1, 2, 3, 4, 4, 4], 2, None),
    ([0, 1, 0, 1, 0, 1], [0, 0, 0, 1, 2, 3, 4, 4, 4], [2], None),
    ([0, 1, 0, 1, 0, 1], [[0, 0, 0, 1, 2, 3, 4, 4, 4]], [2], None),
    ([0, 1, 0, 1, 0, 1], [0, 0, 0, 1, 2, 3, 4, 4, 4], [2, 1], ValueError),
    ([0, 1, 0, 1, 0], [0, 0, 0, 1, 2, 3, 4, 4, 4], 2, Exception),  # Nonconforming number of knots and control points
    ([0, 1, 0, 1, 0, 1], [[0, 0, 0, 1, 2, 3, 4, 4, 4]], 2, None),
    ([0] * (2*3), [[0, 0, 1, 1], [0, 0, 0, 1, 1, 1]], [1, 2], None),
    ([1] * (2*3), [[0, 0, 1, 1], [0, 0, 0, 1, 1, 1]], [1, 2], None),
    ([[cp] for cp in [1] * (2*3)], [[0, 0, 1, 1], [0, 0, 0, 1, 1, 1]], [1, 2], None),
    ([[0, 1] for i in range(2*3)], [[0, 0, 1, 1], [0, 0, 0, 1, 1, 1]], [1, 2], None),
    ([[0, 1, 2] for i in range(2*3)], [[0, 0, 1, 1], [0, 0, 0, 1, 1, 1]], [1, 2], None),
]


@pytest.mark.parametrize("control_points, knot_vectors, degrees, exception", test_param)
def test_bspline_construction(control_points, knot_vectors, degrees, exception):
    # B-spline built from parameters: coefficients, knot vectors and degrees
    if not exception:
        splinterpy.BSpline.init_from_param(control_points, knot_vectors, degrees)
    else:
        with pytest.raises(exception) as e:
            splinterpy.BSpline.init_from_param(control_points, knot_vectors, degrees)
