/*
This file is part of the Multivariate Splines library.
Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#ifndef DATATABLE_H
#define DATATABLE_H

#include <set>
#include "include/datasample.h"

namespace MultivariateSplines
{

// TODO: namespace for interpolation classes

#define SAVE_DOUBLE_PRECISION 17

/*
 * DataTable is a class for storing multidimensional data samples (x,y).
 * The samples are stored in a continuously sorted table.
 */
class DataTable
{
public:
    DataTable();
    DataTable(bool allowDuplicates);
    DataTable(bool allowDuplicates, bool allowIncompleteGrid);

    /*
     * Functions for adding a sample (x,y)
     */
    void addSample(const DataSample &sample);
    void addSample(double x, double y);
    void addSample(std::vector<double> x, double y);
    void addSample(DenseVector x, double y);

    std::multiset<DataSample>::const_iterator cbegin() const;
    std::multiset<DataSample>::const_iterator cend() const;

    unsigned int getNumVariables() const {return numVariables;}
    unsigned int getNumSamples() const {return samples.size();}

    std::vector<std::set<double>> getGrid() const { return grid; }

    /*
     * Backwards compatibility
     */
    std::vector< std::vector<double> > getTableX() const;
    std::vector< std::vector<double> > getTableY() const;
    std::vector<double> getVectorY() const;

    std::vector< std::vector<double> > transposeTable(const std::vector< std::vector<double> > &table) const;
    std::vector< std::vector<double> > getTransposedTableX() const;
    std::vector< std::vector<double> > getTransposedTableY() const;

    /*
     * Save and load functionality
     */
    void save(std::string fileName) const;  // Throws std::ios_base::failure
    void load(std::string fileName);        // Throws std::ios_base::failure

    /*
     * Debug
     */
    void printSamples() const;  // Print point and values for all samples
    void printGrid() const;     // Print the grid (that we know of so far)

    bool isGridComplete() const;

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
};

} // namespace MultivariateSplines

#endif // DATATABLE_H
