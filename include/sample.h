/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_SAMPLE_H
#define SPLINTER_SAMPLE_H

#include <set>
#include "samplepoint.h"

#include <ostream>

namespace SPLINTER
{

/*
 * Sample is a class for storing multidimensional observations (x,y).
 * The samples are stored in a continuously sorted table.
 */
class SPLINTER_API Sample
{
public:
    Sample();
    Sample(bool allowDuplicates);
    Sample(bool allowDuplicates, bool allowIncompleteGrid);
    Sample(const char *fileName);
    Sample(const std::string fileName); // Load Sample from file

    /*
     * Functions for adding a sample point (x,y)
     */
    void addSamplePoint(const SamplePoint &sample);
    void addSamplePoint(double x, double y);
    void addSamplePoint(std::vector<double> x, double y);
    void addSamplePoint(DenseVector x, double y);

    /*
     * Getters
     */
    std::multiset<SamplePoint>::const_iterator cbegin() const;
    std::multiset<SamplePoint>::const_iterator cend() const;

    unsigned int getNumVariables() const {return numVariables;}
    unsigned int size() const {return samples.size();}
    const std::multiset<SamplePoint>& getSamples() const {return samples;}

    std::vector<std::set<double>> getGrid() const { return grid; }
    std::vector< std::vector<double> > getTableX() const;
    std::vector<double> getVectorY() const;
    
    bool isGridComplete() const;

    void save(const std::string fileName) const;

private:
    bool allowDuplicates;
    bool allowIncompleteGrid;
    unsigned int numDuplicates;
    unsigned int numVariables;

    std::multiset<SamplePoint> samples;
    std::vector< std::set<double> > grid;

    void initDataStructures(); // Initialise grid to be a std::vector of xDim std::sets
    unsigned int getNumGridPointsRequired() const;

    void recordGridPoint(const SamplePoint &sample);

    // Used by functions that require the grid to be complete before they start their operation
    // This function prints a message and exits the program if the grid is not complete.
    void gridCompleteGuard() const;

    void load(const std::string fileName);

    friend class Serializer;
    friend bool operator==(const Sample &lhs, const Sample &rhs);
};

Sample operator+(const Sample &lhs, const Sample &rhs);
Sample operator-(const Sample &lhs, const Sample &rhs);

} // namespace SPLINTER

#endif // SPLINTER_SAMPLE_H
