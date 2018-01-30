/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "json_parser.h"
#include "bspline.h"
#include "utilities.h"
#include "fstream"


namespace SPLINTER {

void save_to_json(const BSpline &bspline, const std::string &filename)
{
    std::ofstream ofs(filename);
    nlohmann::json json;

    auto dim_x = bspline.get_dim_x();
    auto dim_y = bspline.get_dim_y();
    auto control_points = eig_mat_to_std_vec_vec(bspline.get_control_points().transpose());
    auto knot_vectors = bspline.get_knot_vectors();
    auto degrees = bspline.get_basis_degrees();

    json["dim_x"] = dim_x;
    json["dim_y"] = dim_y;
    json["degrees"] = degrees;

    for (unsigned int i = 0; i < dim_x; ++i) {
        std::string kv_str = "knot_vector_" + std::to_string(i);
        json[kv_str] = knot_vectors.at(i);
    }

    for (unsigned int i = 0; i < dim_y; ++i) {
        std::string cp_str = "control_points_" + std::to_string(i);
        json[cp_str] = control_points.at(i);
    }

    ofs << json;
}

BSpline load_from_json(const std::string &filename)
{
    std::ifstream ifs(filename);
    nlohmann::json json;
    ifs >> json;

    unsigned int dim_x = json["dim_x"];
    unsigned int dim_y = json["dim_y"];

    std::vector<std::vector<double>> control_points_as_read;
    std::vector<std::vector<double>> control_points;
    std::vector<std::vector<double>> knot_vectors;
    std::vector<unsigned int> degrees;

    auto deg = json["degrees"];
    for (auto &it : deg)
        degrees.push_back(it);

    for (unsigned int i = 0; i < dim_y; ++i) {
        std::vector<double> cp;
        std::string e = "control_points_" + std::to_string(i);
        auto c = json[e];
        for (auto &it : c)
            cp.push_back(it);
        control_points_as_read.push_back(cp);
    }

    // Create list of dim_y-dimensional control points
    for (unsigned int i = 0; i < control_points_as_read.at(0).size(); ++i) {
        std::vector<double> cp;
        for (unsigned int j = 0; j < dim_y; ++j) {
            auto cp_ij = control_points_as_read.at(j).at(i);
            cp.push_back(cp_ij);
        }
        control_points.push_back(cp);
    }

    for (unsigned int i = 0; i < dim_x; ++i) {
        std::vector<double> knots;
        std::string e = "knot_vector_" + std::to_string(i);
        auto k = json[e];
        for (auto &it : k)
            knots.push_back(it);
        knot_vectors.push_back(knots);
    }

    return BSpline(control_points, knot_vectors, degrees);
}

} // namespace SPLINTER