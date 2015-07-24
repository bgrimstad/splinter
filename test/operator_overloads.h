/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_OPERATOR_OVERLOADS_H
#define SPLINTER_OPERATOR_OVERLOADS_H

#include <vector>
#include <set>
#include <datasample.h>
#include <datatable.h>
#include <bspline.h>
#include <radialbasisfunction.h>
#include <polynomialregression.h>
#include <pspline.h>

namespace SPLINTER
{

/*
 * Comparison operators
 */
bool operator==(const DataTable &lhs, const DataTable &rhs);
bool operator==(const DataSample &lhs, const DataSample &rhs);
bool operator==(const BSpline &lhs, const BSpline &rhs);
bool operator==(const PSpline &lhs, const PSpline &rhs);
bool operator==(const RadialBasisFunction &lhs, const RadialBasisFunction &rhs);
bool operator==(const PolynomialRegression &lhs, const PolynomialRegression &rhs);
template <class T>
bool operator==(const std::vector<T> &lhs, const std::vector<T> &rhs);
template <class T>
bool operator==(const std::set<T> &lhs, const std::set<T> &rhs);
template <class T>
bool operator==(const std::multiset<T> &lhs, const std::multiset<T> &rhs);

template <class T>
bool operator!=(const T &lhs, const T &rhs);

/*
 * Output stream operator
 */
std::ostream &operator<<(std::ostream &out, const DataSample &sample);
std::ostream &operator<<(std::ostream &out, const DataTable &table);
template <class T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &obj);
template <class T>
std::ostream &operator<<(std::ostream &out, const std::set<T> &obj);
template <class T>
std::ostream &operator<<(std::ostream &out, const std::multiset<T> &obj);

}

#endif // SPLINTER_OPERATOR_OVERLOADS_H
