import pytest
from os import sys, path, remove
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import splinterpy

# Only for dev purposes
import os.path
if os.path.isdir("/home/bjarne/"):
    splinterpy.load("/home/bjarne/Code/splinter/build/debug/libsplinter-4-0.so")
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
        splinterpy.BSpline.from_param(degrees, knot_vectors, control_points)
    else:
        with pytest.raises(exception) as e:
            splinterpy.BSpline.from_param(degrees, knot_vectors, control_points)


@pytest.mark.parametrize("control_points, knot_vectors, degrees, exception", test_param)
def test_bspline_getters(control_points, knot_vectors, degrees, exception):

    if exception:
        return

    # B-spline built from parameters: coefficients, knot vectors and degrees
    bspline = splinterpy.BSpline.from_param(degrees, knot_vectors, control_points)

    # Compare control points
    control_points_get = bspline.get_control_points()
    for i in range(len(control_points)):
        cp = control_points[i]
        cp_get = control_points_get[i]
        if not isinstance(cp, list):
            assert([cp] == cp_get)
        else:
            assert(cp == cp_get)

    # Compare knot vectors
    knot_vectors_get = bspline.get_knot_vectors()
    if not any(isinstance(kv, list) for kv in knot_vectors):
        assert(knot_vectors_get == [knot_vectors])
    else:
        assert(knot_vectors_get == knot_vectors)

    # Compare degrees
    degrees_get = bspline.get_degrees()
    if not isinstance(degrees, list):
        assert(degrees_get == [degrees])
    else:
        assert(degrees_get == degrees)


def test_bspline_save_load():
    control_points = [0, 1, 0, 1, 0]
    knots = [0, 0, 1, 2, 3, 4, 4]
    degree = 1
    bs1 = splinterpy.BSpline.from_param(degree, knots, control_points)

    filename = "bspline.json"

    try:
        bs1.to_json(filename)

        bs2 = splinterpy.BSpline.from_json(filename)

        assert(all(x == y for x, y in zip(bs1.get_degrees(), bs2.get_degrees())))

    finally:
        remove(filename)
