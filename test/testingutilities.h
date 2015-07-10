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

bool compareFunctions(const Function &exact, const Function &approx, const std::vector<std::vector<double>> &points);

bool compareFunctions(const Function &exact, const Function &approx, const std::vector<std::vector<double>> &points, double val_epsilon, double one_norm_epsilon, double two_norm_epsilon, double inf_norm_epsilon);

bool compareFunctions(const Function &exact, const Function &approx);

bool compareBSplines(const BSpline &left, const BSpline &right);

bool compareDataTables(DataTable &a, DataTable &b);

DataTable sample(const Function &func, std::vector<std::vector<double>> &points);

std::vector<double> linspace(double start, double stop, unsigned int points);

std::vector<std::vector<double>> linspace(std::vector<double> start, std::vector<double> end, std::vector<unsigned int> points);

std::vector<std::vector<double>> linspace(int dim, double start, double end, unsigned int points);

// Returns a default linspace of dim dim
std::vector<std::vector<double>> linspace(int dim);

std::vector<std::vector<double>> linspace(int dim, unsigned int pointsPerDim);

double sixHumpCamelBack(DenseVector x);

double oneNorm(const DenseMatrix &m);

double twoNorm(const DenseMatrix &m);

double maxNorm(const DenseMatrix &m);

// returns log(x) in base base
double log(double base, double x);

std::vector<double> denseToVec(const DenseVector &dense);

DenseVector vecToDense(const std::vector<double> &vec);

std::string pretty_print(const DenseVector &denseVec);

} // namespace SPLINTER

#endif // SPLINTER_TESTINGUTILITIES_H