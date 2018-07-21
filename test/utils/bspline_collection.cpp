/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline_collection.h"
#include "utilities.h"
#include <utils/test_utils.h>


namespace SPLINTER {

std::vector<std::vector<double>> get_random_control_points(unsigned int num_points, unsigned int dim_y)
{
    std::vector<std::vector<double>> control_points;
    for (unsigned int i = 0; i < dim_y; ++i) {
        control_points.push_back(get_random_vector(num_points, 1234+i));
    }
    return transpose_vec_vec(control_points);
}

unsigned int compute_num_control_points(std::vector<std::vector<double>> knot_vectors, std::vector<unsigned int> degrees)
{
    unsigned int num_cp = 1;

    if (knot_vectors.size() != degrees.size()) {
        throw Exception("compute_num_control_points:: Incompatible vector sizes!");
    }

    for (size_t i = 0; i < knot_vectors.size(); ++i) {
        auto num_knots = knot_vectors.at(i).size();
        auto num_basis_functions = num_knots - degrees.at(i) - 1;
        if (num_basis_functions < 1) {
            throw Exception("compute_num_control_points:: Insufficient number of knots!");
        }
        num_cp *= num_basis_functions;
    }
    return num_cp;
}

std::vector<BSpline> get_bspline_collection()
{
    std::vector<BSpline> collection;

    std::vector<double> eqdist_knots = {0, 1, 2, 3, 4, 5};

    // f : R^dim_x -> R^dim_y of various degrees
    {
        for (unsigned int dim_x = 1; dim_x < 4; dim_x++) {
            for (unsigned int dim_y = 1; dim_y < 4; dim_y++) {
                for (unsigned int deg = 1; deg < 4; deg++) {
                    std::vector<unsigned int> degrees;
                    std::vector<std::vector<double>> knot_vectors;
                    for (unsigned int i = 0; i < dim_x; i++) {
                        degrees.push_back(deg);
                        knot_vectors.push_back(eqdist_knots);
                    }
                    auto num_control_points = compute_num_control_points(knot_vectors, degrees);
                    auto control_points = get_random_control_points(num_control_points, dim_y);

                    BSpline bspline(degrees, knot_vectors, control_points);
                    collection.push_back(bspline);
                }
            }
        }
    }

    return collection;
}

} // namespace SPLINTER