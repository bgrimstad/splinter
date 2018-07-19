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
#include <initializer_list>

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

void DataTable::addSample(double x, double y)
{
    addSample(DataPoint(x, y));
}

void DataTable::addSample(const std::vector<double> &x, double y)
{
    addSample(DataPoint(x, y));
}

void DataTable::addSample(double x, const std::vector<double> &y)
{
    addSample(DataPoint(x, y));
}

void DataTable::addSample(const std::vector<double> &x, const std::vector<double> &y)
{
    addSample(DataPoint(x, y));
}

void DataTable::addSample(const DataPoint &sample)
{
    if (getNumSamples() == 0)
    {
        dim_x = sample.getDimX();
        dim_y = sample.getDimY();
        initDataStructures();
    }

    if (sample.getDimX() != dim_x || sample.getDimY() != dim_y) {
        throw Exception("Datatable::addSample: Dimension of new sample is inconsistent with previous samples!");
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

    recordGridPoint(sample);
}

void DataTable::addSample(std::initializer_list<DataPoint> samples)
{
	for (auto& sample : samples)
	{
		addSample(sample);
	}
}

void DataTable::recordGridPoint(const DataPoint &sample)
{
    for (unsigned int i = 0; i < getDimX(); i++)
    {
        grid.at(i).insert(sample.getX().at(i));
    }
}

unsigned int DataTable::getNumSamplesRequired() const
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

bool DataTable::isGridComplete() const
{
    return samples.size() > 0 && samples.size() - num_duplicates == getNumSamplesRequired();
}

void DataTable::initDataStructures()
{
    for (unsigned int i = 0; i < getDimX(); i++)
    {
        grid.push_back(std::set<double>());
    }
}

void DataTable::gridCompleteGuard() const
{
    if (!(isGridComplete() || allow_incomplete_grid))
    {
        throw Exception("DataTable::gridCompleteGuard: The grid is not complete yet!");
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
std::vector<std::vector<double>> DataTable::getTableX() const
{
    gridCompleteGuard();

    std::vector<std::vector<double>> table(dim_x, std::vector<double>(getNumSamples()));

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
std::vector<std::vector<double>> DataTable::getTableY() const
{
    std::vector<std::vector<double>> table(dim_y, std::vector<double>(getNumSamples()));
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
    if(lhs.getDimX() != rhs.getDimX()) {
        throw Exception("operator+(DataTable, DataTable): trying to add two DataTable's of different dimensions!");
    }

    DataTable result;
    for(auto it = lhs.cbegin(); it != lhs.cend(); it++) {
        result.addSample(*it);
    }
    for(auto it = rhs.cbegin(); it != rhs.cend(); it++) {
        result.addSample(*it);
    }

    return result;
}

DataTable operator-(const DataTable &lhs, const DataTable &rhs)
{
    if(lhs.getDimX() != rhs.getDimX()) {
        throw Exception("operator-(DataTable, DataTable): trying to subtract two DataTable's of different dimensions!");
    }

    DataTable result;
    auto rhsSamples = rhs.getSamples();
    // Add all samples from lhs that are not in rhs
    for(auto it = lhs.cbegin(); it != lhs.cend(); it++) {
        if(rhsSamples.count(*it) == 0) {
            result.addSample(*it);
        }
    }

    return result;
}

} // namespace SPLINTER
