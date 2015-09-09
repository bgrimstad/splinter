# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


__all__ = ["splinter", "datatable", "bspline", "pspline", "radialbasisfunction", "polynomialregression"]

from .datatable import DataTable
from .bspline import BSpline
from .pspline import PSpline
from .radialbasisfunction import RadialBasisFunction, RBFType
from .polynomialregression import PolynomialRegression

splinter.DataTable = datatable.DataTable
splinter.BSpline = bspline.BSpline
splinter.PSpline = pspline.PSpline
splinter.RadialBasisFunction = radialbasisfunction.RadialBasisFunction
splinter.RBFType = radialbasisfunction.RBFType
splinter.PolynomialRegression = polynomialregression.PolynomialRegression

from .splinter import *
try:
	load()
except Exception as e:
	print(e)