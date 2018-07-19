/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "data_table2.h"
#include "json_parser.h"
#include <set>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <limits>


namespace SPLINTER
{

DataTable2::DataTable2()
    : _dim_x(1),
      _dim_y(1)
{
}

// TODO: add return value (iterator)
void DataTable2::add_sample(const DataPoint &sample)
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

/*
 * Get table of samples x-values: i.e. table[i][j] is the value of input i at sample j
 */
std::vector<std::vector<double>> DataTable2::get_table_x() const
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

/*
 * Get table of sampled y-values: i.e. table[i][j] is the value of output i at sample j
 */
std::vector<std::vector<double>> DataTable2::get_table_y() const
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

} // namespace SPLINTER
