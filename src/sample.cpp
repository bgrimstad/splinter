/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "sample.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <serializer.h>

namespace SPLINTER
{

Sample::Sample()
    : Sample(false, false)
{
}

Sample::Sample(bool allowDuplicates)
    : Sample(allowDuplicates, false)
{
}

Sample::Sample(bool allowDuplicates, bool allowIncompleteGrid)
    : allowDuplicates(allowDuplicates),
      allowIncompleteGrid(allowIncompleteGrid),
      numDuplicates(0),
      numVariables(0)
{
}

Sample::Sample(const char *fileName)
    : Sample(std::string(fileName))
{
}

Sample::Sample(const std::string fileName)
{
    load(fileName);
}

void Sample::addSample(double x, double y)
{
    addSample(SamplePoint(x, y));
}

void Sample::addSample(std::vector<double> x, double y)
{
    addSample(SamplePoint(x, y));
}

void Sample::addSample(DenseVector x, double y)
{
    addSample(SamplePoint(x, y));
}

void Sample::addSample(const SamplePoint &sample)
{
    if (getNumSamples() == 0)
    {
        numVariables = sample.getDimX();
        initDataStructures();
    }

    if(sample.getDimX() != numVariables) {
        throw Exception("Sample::addSample: Dimension of new sample is inconsistent with previous samples!");
    }

    // Check if the sample has been added already
    if (samples.count(sample) > 0)
    {
        if (!allowDuplicates)
        {
#ifndef NDEBUG
            std::cout << "Discarding duplicate sample because allowDuplicates is false!" << std::endl;
            std::cout << "Initialise with Sample(true) to set it to true." << std::endl;
#endif // NDEBUG

            return;
        }

        numDuplicates++;
    }

    samples.insert(sample);

    recordGridPoint(sample);
}

void Sample::recordGridPoint(const SamplePoint &sample)
{
    for (unsigned int i = 0; i < getNumVariables(); i++)
    {
        grid.at(i).insert(sample.getX().at(i));
    }
}

unsigned int Sample::getNumSamplesRequired() const
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

bool Sample::isGridComplete() const
{
    return samples.size() > 0 && samples.size() - numDuplicates == getNumSamplesRequired();
}

void Sample::initDataStructures()
{
    for (unsigned int i = 0; i < getNumVariables(); i++)
    {
        grid.push_back(std::set<double>());
    }
}

void Sample::gridCompleteGuard() const
{
    if (!(isGridComplete() || allowIncompleteGrid))
    {
        throw Exception("Sample::gridCompleteGuard: The grid is not complete yet!");
    }
}

void Sample::save(const std::string fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void Sample::load(const std::string fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

/*
 * Getters for iterators
 */
std::multiset<SamplePoint>::const_iterator Sample::cbegin() const
{
    return samples.cbegin();
}

std::multiset<SamplePoint>::const_iterator Sample::cend() const
{
    return samples.cend();
}

/*
 * Get table of samples x-values,
 * i.e. table[i][j] is the value of variable i at sample j
 */
std::vector< std::vector<double> > Sample::getTableX() const
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

        for (unsigned int j = 0; j < numVariables; j++)
        {
            table.at(j).at(i) = x.at(j);
        }
        i++;
    }

    return table;
}

// Get vector of y-values
std::vector<double> Sample::getVectorY() const
{
    std::vector<double> y;
    for (std::multiset<SamplePoint>::const_iterator it = cbegin(); it != cend(); ++it)
    {
        y.push_back(it->getY());
    }
    return y;
}

Sample operator+(const Sample &lhs, const Sample &rhs)
{
    if(lhs.getNumVariables() != rhs.getNumVariables()) {
        throw Exception("operator+(Sample, Sample): trying to add two Sample's of different dimensions!");
    }

    Sample result;
    for(auto it = lhs.cbegin(); it != lhs.cend(); it++) {
        result.addSample(*it);
    }
    for(auto it = rhs.cbegin(); it != rhs.cend(); it++) {
        result.addSample(*it);
    }

    return result;
}

Sample operator-(const Sample &lhs, const Sample &rhs)
{
    if(lhs.getNumVariables() != rhs.getNumVariables()) {
        throw Exception("operator-(Sample, Sample): trying to subtract two Sample's of different dimensions!");
    }

    Sample result;
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
