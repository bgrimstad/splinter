/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "data_table.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <serializer.h>

namespace SPLINTER
{

DataTable::DataTable()
    : DataTable(false, false)
{
}

DataTable::DataTable(bool allowDuplicates)
    : DataTable(allowDuplicates, false)
{
}

DataTable::DataTable(bool allowDuplicates, bool allowIncompleteGrid)
    : allow_duplicates(allowDuplicates),
      allow_incomplete_grid(allowIncompleteGrid),
      num_duplicates(0),
      dim_x(0),
      dim_y(0)
{
}

DataTable::DataTable(const char *fileName)
    : DataTable(std::string(fileName))
{
}

DataTable::DataTable(const std::string &fileName)
{
    load(fileName);
}

void DataTable::add_sample(double x, double y)
{
    add_sample(DataPoint(x, y));
}

void DataTable::add_sample(const std::vector<double> &x, double y)
{
    add_sample(DataPoint(x, y));
}

void DataTable::add_sample(double x, const std::vector<double> &y)
{
    add_sample(DataPoint(x, y));
}

void DataTable::add_sample(const std::vector<double> &x, const std::vector<double> &y)
{
    add_sample(DataPoint(x, y));
}

void DataTable::add_sample(const DataPoint &sample)
{
    if (get_num_samples() == 0)
    {
        dim_x = sample.get_dim_x();
        dim_y = sample.get_dim_y();
        init_data_structures();
    }

    if (sample.get_dim_x() != dim_x || sample.get_dim_y() != dim_y) {
        throw Exception("Datatable::add_sample: Dimension of new sample is inconsistent with previous samples!");
    }

    // Check if the sample has been added already
    if (samples.count(sample) > 0)
    {
        if (!allow_duplicates)
        {
#ifndef NDEBUG
            std::cout << "Discarding duplicate sample because allow_duplicates is false!" << std::endl;
            std::cout << "Initialise with DataTable(true) to set it to true." << std::endl;
#endif // NDEBUG

            return;
        }

        num_duplicates++;
    }

    samples.insert(sample);

    record_grid_point(sample);
}

void DataTable::record_grid_point(const DataPoint &sample)
{
    for (unsigned int i = 0; i < get_dim_x(); i++)
    {
        grid.at(i).insert(sample.getX().at(i));
    }
}

unsigned int DataTable::get_num_samples_required() const
{
    unsigned long samplesRequired = 1;
    unsigned int i = 0;
    for (auto &variable : grid)
    {
        samplesRequired *= variable.size();
        i++;
    }

    return (i > 0 ? samplesRequired : (unsigned long) 0);
}

bool DataTable::is_grid_complete() const
{
    return samples.size() > 0 && samples.size() - num_duplicates == get_num_samples_required();
}

void DataTable::init_data_structures()
{
    for (unsigned int i = 0; i < get_dim_x(); i++)
    {
        grid.push_back(std::set<double>());
    }
}

void DataTable::grid_complete_guard() const
{
    if (!(is_grid_complete() || allow_incomplete_grid))
    {
        throw Exception("DataTable::grid_complete_guard: The grid is not complete yet!");
    }
}

void DataTable::save(const std::string &fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void DataTable::load(const std::string &fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

/*
 * Getters for iterators
 */
std::multiset<DataPoint>::const_iterator DataTable::cbegin() const
{
    return samples.cbegin();
}

std::multiset<DataPoint>::const_iterator DataTable::cend() const
{
    return samples.cend();
}

/*
 * Get table of samples x-values: i.e. table[i][j] is the value of input i at sample j
 */
std::vector<std::vector<double>> DataTable::get_table_x() const
{
    grid_complete_guard();

    std::vector<std::vector<double>> table(dim_x, std::vector<double>(get_num_samples()));

    unsigned int i = 0;
    for (auto &sample : samples) {
        auto x = sample.getX();
        for (unsigned int j = 0; j < dim_x; j++)
            table.at(j).at(i) = x.at(j);
        i++;
    }

    return table;
}

/*
 * Get table of sampled y-values: i.e. table[i][j] is the value of output i at sample j
 */
std::vector<std::vector<double>> DataTable::get_table_y() const
{
    std::vector<std::vector<double>> table(dim_y, std::vector<double>(get_num_samples()));
    unsigned int i = 0;
    for (const auto &sample : samples) {
        auto y = sample.getY();
        for (unsigned int j = 0; j < dim_y; ++j)
            table.at(j).at(i) = y.at(j);
        i++;
    }

    return table;
}

DataTable operator+(const DataTable &lhs, const DataTable &rhs)
{
    if(lhs.get_dim_x() != rhs.get_dim_x()) {
        throw Exception("operator+(DataTable, DataTable): trying to add two DataTable's of different dimensions!");
    }

    DataTable result;
    for(auto it = lhs.cbegin(); it != lhs.cend(); it++) {
        result.add_sample(*it);
    }
    for(auto it = rhs.cbegin(); it != rhs.cend(); it++) {
        result.add_sample(*it);
    }

    return result;
}

DataTable operator-(const DataTable &lhs, const DataTable &rhs)
{
    if(lhs.get_dim_x() != rhs.get_dim_x()) {
        throw Exception("operator-(DataTable, DataTable): trying to subtract two DataTable's of different dimensions!");
    }

    DataTable result;
    auto rhsSamples = rhs.get_samples();
    // Add all samples from lhs that are not in rhs
    for(auto it = lhs.cbegin(); it != lhs.cend(); it++) {
        if(rhsSamples.count(*it) == 0) {
            result.add_sample(*it);
        }
    }

    return result;
}

} // namespace SPLINTER
