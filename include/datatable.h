/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_DATATABLE_H
#define SPLINTER_DATATABLE_H

#include <set>
#include "datasample.h"

#include <ostream>

namespace SPLINTER
{

/*
 * DataTable is a class for storing multidimensional data samples (x,y).
 * The samples are stored in a continuously sorted table.
 */
class API DataTable
{
public:
    DataTable();
    DataTable(bool allowDuplicates);
    DataTable(bool allowDuplicates, bool allowIncompleteGrid);
    DataTable(const char *fileName);
    DataTable(const std::string fileName); // Load DataTable from file

    /*
     * Functions for adding a sample (x,y)
     */
    void addSample(const DataSample &sample);
    void addSample(double x, double y);
    void addSample(std::vector<double> x, double y);
    void addSample(DenseVector x, double y);

    /*
     * Getters
     */
    std::multiset<DataSample>::const_iterator cbegin() const;
    std::multiset<DataSample>::const_iterator cend() const;

    unsigned int getNumVariables() const {return numVariables;}
    unsigned int getNumSamples() const {return samples.size();}
    const std::multiset<DataSample>& getSamples() const {return samples;}

    std::vector<std::set<double>> getGrid() const { return grid; }
    std::vector< std::vector<double> > getTableX() const;
    std::vector<double> getVectorY() const;

    /*
     * Debug
     */
    template <typename _Char>
    void printSamples(std::basic_ostream<_Char>& out) const
    {
        for (auto &sample : samples)
        {
            out << sample << std::endl;
        }
    }

    template <typename _Char>
    void printGrid(std::basic_ostream<_Char>& out) const
    {
        out << "===== Printing grid =====" << std::endl;

        unsigned int i = 0;
        for (auto &variable : grid)
        {
            out << 'x' << i++ << '(' << variable.size() << "): ";
            for (double value : variable)
            {
                out << value << ' ';
            }
            out << std::endl;
        }

        out << "Unique samples added: " << samples.size() << std::endl;
        out << "Samples required: " << getNumSamplesRequired() << std::endl;
    }

    bool isGridComplete() const;

    void save(const std::string fileName) const;

private:
    bool allowDuplicates;
    bool allowIncompleteGrid;
    unsigned int numDuplicates;
    unsigned int numVariables;

    std::multiset<DataSample> samples;
    std::vector< std::set<double> > grid;

    void initDataStructures(); // Initialise grid to be a std::vector of xDim std::sets
    unsigned int getNumSamplesRequired() const;

    void recordGridPoint(const DataSample &sample);

    // Used by functions that require the grid to be complete before they start their operation
    // This function prints a message and exits the program if the grid is not complete.
    void gridCompleteGuard() const;

    void load(const std::string fileName);

    friend class Serializer;
};

} // namespace SPLINTER

#endif // SPLINTER_DATATABLE_H