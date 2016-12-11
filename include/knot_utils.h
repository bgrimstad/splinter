/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_KNOT_UTILS_H
#define SPLINTER_KNOT_UTILS_H

#include <vector>

namespace SPLINTER
{

// Knot vector construction
std::vector<double> knotVectorEquidistantNotClamped(const std::vector<double> &values,
                                                    unsigned int degree,
                                                    unsigned int numBasisFunctions);

std::vector<double> knotVectorMovingAverage(const std::vector<double> &values,
                                            unsigned int degree);

std::vector<double> knotVectorEquidistant(const std::vector<double> &values,
                                          unsigned int degree,
                                          unsigned int numBasisFunctions);

} // namespace SPLINTER

#endif // SPLINTER_KNOT_UTILS_H
