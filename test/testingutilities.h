/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef MS_TESTINGUTILITIES_H
#define MS_TESTINGUTILITIES_H

#include <datatable.h>
#include <bspline.h>
#include <generaldefinitions.h>

namespace MultivariateSplines
{

bool equalsWithinRange(double a, double b, double margin = 0.0);

bool compareBSplines(BSpline &bs, const BSpline &bs_orig);

bool compareDataTables(DataTable &a, DataTable &b);

std::vector<double> linspace(double start, double stop, unsigned int points);

double sixHumpCamelBack(DenseVector x);

} // namespace MultivariateSplines

#endif // MS_TESTINGUTILITIES_H
