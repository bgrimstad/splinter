/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "data_table_collection.h"
#include "utilities.h"
#include "test_utils.h"
#include "test_function_utils.h"


namespace SPLINTER {

std::vector<DataTable> get_data_table_collection()
{
    std::vector<DataTable> collection;

    double x_start = -10;
    double x_stop = 10;
    unsigned int num_points = 10;  // Number of points per input dimension

    // Create data tables with random samples of various dimensions
    auto counter = 1;
    for (unsigned int dim_x = 1; dim_x <= 3; dim_x++) {
        for (unsigned int dim_y = 1; dim_y <= 3; dim_y++) {
            // New DataTable
            auto table = DataTable();

            // Create grid to sample on
            auto x_start_vector = std::vector<double>(dim_x, x_start);
            auto x_stop_vector = std::vector<double>(dim_x, x_stop);
            auto num_points_vector = std::vector<unsigned int>(dim_x, num_points);
            auto points = multi_linspace(x_start_vector, x_stop_vector, num_points_vector);

            for (auto &x : points) {
                // Create random sample values
                auto seed = dim_x*dim_y*123*counter;
                ++counter;
                auto y = get_random_vector(dim_y, seed);
                table.add_sample(x, y);
            }

            collection.emplace_back(table);
        }
    }

    return collection;
}

} // namespace SPLINTER