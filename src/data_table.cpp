/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "data_table.h"
#include "json_parser.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <initializer_list>


namespace SPLINTER
{

DataTable::DataTable()
    : DataTable(false, false)
{
}

DataTable::DataTable(bool allow_duplicates)
    : DataTable(allow_duplicates, false)
{
}

DataTable::DataTable(bool allow_duplicates, bool allow_incomplete_grid)
    : _allow_duplicates(allow_duplicates),
      _allow_incomplete_grid(allow_incomplete_grid),
      _num_duplicates(0),
      _dim_x(0),
      _dim_y(0)
{
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
        _dim_x = sample.get_dim_x();
        _dim_y = sample.get_dim_y();
        init_data_structures();
    }

    if (sample.get_dim_x() != _dim_x || sample.get_dim_y() != _dim_y) {
        throw Exception("DataTable::add_sample: Dimension of new sample is inconsistent with previous samples!");
    }

    // Check if the sample has been added already
    if (samples.count(sample) > 0)
    {
        if (!_allow_duplicates)
        {
#ifndef NDEBUG
            std::cout << "Discarding duplicate sample because _allow_duplicates is false!" << std::endl;
            std::cout << "Initialise with DataTable(true) to set it to true." << std::endl;
#endif // NDEBUG

            return;
        }

        _num_duplicates++;
    }

    samples.insert(sample);

    record_grid_point(sample);
}

void DataTable::add_sample(std::initializer_list<DataPoint> samples)
{
	for (auto& sample : samples)
	{
		add_sample(sample);
	}
}

void DataTable::record_grid_point(const DataPoint &sample)
{
    for (unsigned int i = 0; i < get_dim_x(); i++)
    {
        grid.at(i).insert(sample.get_x().at(i));
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
    return samples.size() > 0 && samples.size() - _num_duplicates == get_num_samples_required();
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
    if (!(is_grid_complete() || _allow_incomplete_grid))
    {
        throw Exception("DataTable::grid_complete_guard: The grid is not complete yet!");
    }
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

void DataTable::to_json(const std::string &filename) const {
    SPLINTER::datatable_to_json(*this, filename);
}

DataTable DataTable::from_json(const std::string &filename) {
    return SPLINTER::datatable_from_json(filename);
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
    auto rhs_samples = rhs.get_samples();
    // Add all samples from lhs that are not in rhs
    for(auto it = lhs.cbegin(); it != lhs.cend(); it++) {
        if(rhs_samples.count(*it) == 0) {
            result.add_sample(*it);
        }
    }

    return result;
}

} // namespace SPLINTER
