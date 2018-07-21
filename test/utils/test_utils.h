/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_TEST_UTILS_H
#define SPLINTER_TEST_UTILS_H

#include <data_table.h>
#include <function.h>
#include "definitions.h"
#include <utils/op_overloads.h>
#include <utils/test_function.h>


namespace SPLINTER
{

double rel_error(double exactVal, double approxVal);

std::vector<double> get_random_vector(unsigned int length, unsigned int seed);

// points is a vector where each element is the number of points for that dim
std::vector<std::vector<double>> multi_linspace(std::vector<double> start, std::vector<double> end,
                                                std::vector<unsigned int> points);

// points is the total number of points, not per dim
std::vector<std::vector<double>> multi_linspace(int dim, double start, double end, unsigned int points);

// Returns a default linspace of dim dim
std::vector<std::vector<double>> multi_linspace(int dim);

std::vector<std::vector<double>> multi_linspace(int dim, unsigned int points_per_dim);

double six_hump_camel_back(std::vector<double> x);

/*
 * Matrix norms
 */
double matrix_one_norm(const DenseMatrix &mat);

double matrix_two_norm(const DenseMatrix &mat);

double matrix_inf_norm(const DenseMatrix &mat);

// returns log(x) in base base
double log(double base, double x);

std::string pretty_print(const DenseVector &vec);

} // namespace SPLINTER

#endif // SPLINTER_TEST_UTILS_H
