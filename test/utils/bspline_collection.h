/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_COLLECTION_H
#define SPLINTER_BSPLINE_COLLECTION_H

#include "bspline_builders.h"
#include <vector>

namespace SPLINTER {

std::vector<std::vector<double>> get_random_control_points(unsigned int num_points, unsigned int dim_y);

unsigned int compute_num_control_points(std::vector<std::vector<double>> knot_vectors, std::vector<unsigned int> degrees);

std::vector<BSpline> get_bspline_collection();

} // namespace SPLINTER

#endif //SPLINTER_BSPLINE_COLLECTION_H
