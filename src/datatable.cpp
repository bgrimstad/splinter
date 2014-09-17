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

#include <iomanip> // Set precision used by streams on double values
#include <fstream>

namespace MultivariateSplines
{

/*
 * Workaround for https://gcc.gnu.org/bugzilla/show_bug.cgi?id=52015
 * Basically, std::stod and std::stoi are deactivated on MinGW because of a bug
 */
//#if defined(__MINGW32__) || defined(__MINGW64__)
#include <sstream>
    auto stringToDouble = [](std::string s, std::string::size_type *sz)
    {
        double d;
        std::stringstream ss(s); // Turn the string into a stream
        ss >> d; // Convert
        ss.str("");
        ss << std::setprecision(SAVE_DOUBLE_PRECISION) << d;
        *sz = ss.str().size() + 1;
        return d;
    };
//#else
    //auto stringToDouble = std::stod;
//#endif

DataTable::DataTable()
    : DataTable(false, false)
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
    if(getNumSamples() == 0)
    {
        numVariables = sample.getDimX();
        initDataStructures();
    }

    assert(sample.getDimX() == numVariables); // All points must have the same dimension

    // Check if the sample has been added already
    if(samples.count(sample) > 0)
    {
        if(!allowDuplicates)
        {
            std::cout << "Discarding duplicate sample because allowDuplicates is false!" << std::endl;
            std::cout << "Initialise with DataTable(true) to set it to true." << std::endl;
            return;
        }

        numDuplicates++;
    }

    samples.insert(sample);

    recordGridPoint(sample);
}

void DataTable::recordGridPoint(const DataSample &sample)
{
    for(unsigned int i = 0; i < getNumVariables(); i++)
    {
        grid.at(i).insert(sample.getX().at(i));
    }
}

unsigned int DataTable::getNumSamplesRequired() const
{
    unsigned long samplesRequired = 1;
    unsigned int i = 0;
    for(auto &variable : grid)
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
    for(unsigned int i = 0; i < getNumVariables(); i++)
    {
        grid.push_back(std::set<double>());
    }
}

void DataTable::gridCompleteGuard() const
{
    if(!isGridComplete() && !allowIncompleteGrid)
    {
        std::cout << "The grid is not complete yet!" << std::endl;
        exit(1);
    }
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
    for(unsigned int i = 0; i < numVariables; i++)
    {
        std::vector<double> xi(getNumSamples(), 0.0);
        table.push_back(xi);
    }

    // Fill table with values
    int i = 0;
    for(auto &sample : samples)
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
    for(std::multiset<DataSample>::const_iterator it = cbegin(); it != cend(); ++it)
    {
        y.push_back(it->getY());
    }
    return y;
}

/*****************
 * Save and load *
 *****************/

void DataTable::save(std::string fileName) const
{
    std::ofstream outFile;

    try
    {
        outFile.open(fileName);
    }
    catch(const std::ios_base::failure &e)
    {
        throw;
    }

    // If this function is still alive the file must be open,
    // no need to call is_open()

    // Write header
    outFile << "# Saved DataTable" << '\n';
    outFile << "# Number of samples: " << getNumSamples() << '\n';
    outFile << "# Complete grid: " << (isGridComplete() ? "yes" : "no") << '\n';
    outFile << "# xDim: " << numVariables << '\n';
    outFile << numVariables << " " << 1 << '\n';

    for(auto it = cbegin(); it != cend(); it++)
    {
        for(unsigned int i = 0; i < numVariables; i++)
        {
            outFile << std::setprecision(SAVE_DOUBLE_PRECISION) << it->getX().at(i) << " ";
        }

        outFile << std::setprecision(SAVE_DOUBLE_PRECISION) << it->getY();

        outFile << '\n';
    }

    // close() also flushes
    try
    {
        outFile.close();
    }
    catch(const std::ios_base::failure &e)
    {
        throw;
    }
}

void DataTable::load(std::string fileName)
{
    std::ifstream inFile;

    try
    {
        inFile.open(fileName);
    }
    catch(const std::ios_base::failure &e)
    {
        throw;
    }

    // If this function is still alive the file must be open,
    // no need to call is_open()

    // Skip past comments
    std::string line;

    int nX, nY;
    int state = 0;
    while(std::getline(inFile, line))
    {
        // Look for comment sign
        if(line.at(0) == '#')
            continue;

        // Reading number of dimensions
        if(state == 0)
        {
            std::string::size_type sz = 0;
            nX = (int) round(stringToDouble(line, &sz));
            nY = 1;
            state = 1;
        }

        // Reading samples
        else if(state == 1)
        {
            auto x = std::vector<double>(nX);
            auto y = std::vector<double>(nY);
            std::string::size_type sz = 0;
            for(int i = 0; i < nX; i++)
            {
                line = line.substr(sz);
                x.at(i) = stringToDouble(line, &sz);
            }
            for(int j = 0; j < nY; j++)
            {
                line = line.substr(sz);
                y.at(j) = stringToDouble(line, &sz);
            }

            addSample(x, y.at(0));
        }
    }

    // close() also flushes
    try
    {
        inFile.close();
    }
    catch(const std::ios_base::failure &e)
    {
        throw;
    }
}

/**************
 * Debug code *
 **************/
void DataTable::printSamples() const
{
    for(auto &sample : samples)
    {
        std::cout << sample << std::endl;
    }
}

void DataTable::printGrid() const
{
    std::cout << "===== Printing grid =====" << std::endl;

    unsigned int i = 0;
    for(auto &variable : grid)
    {
        std::cout << "x" << i++ << "(" << variable.size() << "): ";
        for(double value : variable)
        {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Unique samples added: " << samples.size() << std::endl;
    std::cout << "Samples required: " << getNumSamplesRequired() << std::endl;
}

} // namespace MultivariateSplines
