/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "datatable.h"

#include <string>
#include <fstream>
#include <iomanip>

#include <stdexcept>
#include <limits>
#include <serialize.h>

namespace MultivariateSplines
{

DataTable::DataTable()
    : DataTable(false, false)
{
}

DataTable::DataTable(const char *fileName)
    : DataTable(std::string(fileName))
{
}

DataTable::DataTable(bool allowDuplicates)
    : DataTable(allowDuplicates, false)
{
}

DataTable::DataTable(bool allowDuplicates, bool allowIncompleteGrid)
    : allowDuplicates(allowDuplicates),
      allowIncompleteGrid(allowIncompleteGrid),
      numDuplicates(0),
      numVariables(0)
{
}

DataTable::DataTable(const std::string fileName)
{
    load(fileName);
}

void DataTable::addSample(double x, double y)
{
    addSample(DataSample(x, y));
}

void DataTable::addSample(std::vector<double> x, double y)
{
    addSample(DataSample(x, y));
}

void DataTable::addSample(DenseVector x, double y)
{
    addSample(DataSample(x, y));
}

void DataTable::addSample(const DataSample &sample)
{
    if (getNumSamples() == 0)
    {
        numVariables = sample.getDimX();
        initDataStructures();
    }

    assert(sample.getDimX() == numVariables); // All points must have the same dimension

    // Check if the sample has been added already
    if (samples.count(sample) > 0)
    {
        if (!allowDuplicates)
        {
#ifndef NDEBUG
            std::cout << "Discarding duplicate sample because allowDuplicates is false!" << std::endl;
            std::cout << "Initialise with DataTable(true) to set it to true." << std::endl;
#endif // NDEBUG

            return;
        }

        numDuplicates++;
    }

    samples.insert(sample);

    recordGridPoint(sample);
}

void DataTable::recordGridPoint(const DataSample &sample)
{
    for (unsigned int i = 0; i < getNumVariables(); i++)
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
        samplesRequired *= (unsigned long) variable.size();
        i++;
    }

    return (i > 0 ? samplesRequired : (unsigned long) 0);
}

bool DataTable::isGridComplete() const
{
    return samples.size() > 0 && samples.size() - numDuplicates == getNumSamplesRequired();
}

void DataTable::initDataStructures()
{
    for (unsigned int i = 0; i < getNumVariables(); i++)
    {
        grid.push_back(std::set<double>());
    }
}

void DataTable::gridCompleteGuard() const
{
    if (!(isGridComplete() || allowIncompleteGrid))
    {
        throw Exception("DataTable::gridCompleteGuard: The grid is not complete yet!");
    }
}

void DataTable::save(const std::string fileName) const
{
    StreamType stream;
    serialize(allowDuplicates, stream);
    serialize(allowIncompleteGrid, stream);
    serialize(numDuplicates, stream);
    serialize(numVariables, stream);
    serialize(samples, stream);
    serialize(grid, stream);

    save_to_file(fileName, stream);
}

void DataTable::load(const std::string fileName)
{
    StreamType stream = load_from_file(fileName);

    auto it = stream.cbegin();
    allowDuplicates = deserialize<bool>(it, stream.cend());
    allowIncompleteGrid = deserialize<bool>(it, stream.cend());
    numDuplicates = deserialize<unsigned int>(it, stream.cend());
    numVariables = deserialize<unsigned int>(it, stream.cend());
    samples = deserialize<std::multiset<DataSample>>(it, stream.cend());
    grid = deserialize<std::vector<std::set<double>>>(it, stream.cend());
}

/***********
 * Getters *
 ***********/

std::multiset<DataSample>::const_iterator DataTable::cbegin() const
{
    gridCompleteGuard();

    return samples.cbegin();
}

std::multiset<DataSample>::const_iterator DataTable::cend() const
{
    gridCompleteGuard();

    return samples.cend();
}

// Get table of samples x-values,
// i.e. table[i][j] is the value of variable i at sample j
std::vector< std::vector<double> > DataTable::getTableX() const
{
    gridCompleteGuard();

    // Initialize table
    std::vector<std::vector<double>> table;
    for (unsigned int i = 0; i < numVariables; i++)
    {
        std::vector<double> xi(getNumSamples(), 0.0);
        table.push_back(xi);
    }

    // Fill table with values
    int i = 0;
    for (auto &sample : samples)
    {
        std::vector<double> x = sample.getX();

        for(unsigned int j = 0; j < numVariables; j++)
        {
            table.at(j).at(i) = x.at(j);
        }
        i++;
    }

    return table;
}

// Get vector of y-values
std::vector<double> DataTable::getVectorY() const
{
    std::vector<double> y;
    for (std::multiset<DataSample>::const_iterator it = cbegin(); it != cend(); ++it)
    {
        y.push_back(it->getY());
    }
    return y;
}

} // namespace MultivariateSplines
