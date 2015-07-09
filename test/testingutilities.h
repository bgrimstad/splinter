/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_TESTINGUTILITIES_H
#define SPLINTER_TESTINGUTILITIES_H

#include <datatable.h>
#include <function.h>
#include <generaldefinitions.h>
#include <bspline.h>

namespace SPLINTER
{

bool equalsWithinRange(double a, double b, double margin = 0.0);

bool compareFunctions(const Function &f1, const Function &f2);

bool compareBSplines(BSpline &bs, const BSpline &bs_orig);

bool compareDataTables(DataTable &a, DataTable &b);

std::vector<double> linspace(double start, double stop, unsigned int points);

double sixHumpCamelBack(DenseVector x);

double oneNorm(const DenseMatrix &m);

double twoNorm(const DenseMatrix &m);

double maxNorm(const DenseMatrix &m);

// returns log(x) in base base
double log(double base, double x);

} // namespace SPLINTER

#endif // SPLINTER_TESTINGUTILITIES_H