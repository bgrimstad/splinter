/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "data_table.h"
#include <set>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <limits>


namespace SPLINTER
{

// TODO: add return value (iterator)
void DataTable::add_sample(const DataPoint &sample)
{
    if (get_num_samples() == 0) {
        _dim_x = sample.get_dim_x();
        _dim_y = sample.get_dim_y();
    }

    if (sample.get_dim_x() != _dim_x || sample.get_dim_y() != _dim_y) {
        throw Exception("DataTable::add_sample: Dimension of new sample is inconsistent with previous samples!");
    }

    samples.emplace_back(sample);
}

void DataTable::add_sample(std::initializer_list<DataPoint> samples)
{
    for (auto& sample : samples)
    {
        add_sample(sample);
    }
}

/**
 * Get table of samples x-values: i.e. table[i][j] is the value of input i at sample j
 */
std::vector<std::vector<double>> DataTable::get_table_x() const
{
    std::vector<std::vector<double>> table(_dim_x, std::vector<double>(get_num_samples()));

    unsigned int i = 0;
    for (auto &sample : samples) {
        auto x = sample.get_x();
        for (unsigned int j = 0; j < _dim_x; j++)
            table.at(j).at(i) = x.at(j);
        i++;
    }

    return table;
}

/**
 * Get table of sampled y-values: i.e. table[i][j] is the value of output i at sample j
 */
std::vector<std::vector<double>> DataTable::get_table_y() const
{
    std::vector<std::vector<double>> table(_dim_y, std::vector<double>(get_num_samples()));
    unsigned int i = 0;
    for (const auto &sample : samples) {
        auto y = sample.get_y();
        for (unsigned int j = 0; j < _dim_y; ++j)
            table.at(j).at(i) = y.at(j);
        i++;
    }

    return table;
}

/**
 * Check if grid is complete
 */
bool DataTable::is_grid_complete() const
{
    // Return true if there are no samples (ie. the grid is empty)
    if (get_num_samples() == 0)
    {
        return true;
    }

    // Construct grid
    std::vector<std::set<double>> grid;

    for (unsigned int i = 0; i < get_dim_x(); i++)
    {
        grid.emplace_back(std::set<double>());
    }

    for (const auto &sample : samples)
    {
        for (unsigned int i = 0; i < get_dim_x(); i++)
        {
            grid.at(i).insert(sample.get_x().at(i));
        }
    }

    // Compute number of grid points
    unsigned int num_grid_points = 1;
    for (auto &variable : grid)
    {
        num_grid_points *= variable.size();
    }

    // Check that there is a sample at every grid point.
    // First, we count the number of unique sample points using a multiset (uniqueness is based on x-values).
    // If the grid is fully sampled, this number must be equal to the number of grid points.
    std::multiset<DataPoint> unique_samples;
    for (const auto &sample : samples)
    {
        // Check if the sample has been added already
        if (unique_samples.count(sample) == 0)
        {
            unique_samples.insert(sample);
        }
    }

    return unique_samples.size() == num_grid_points;

}

} // namespace SPLINTER
