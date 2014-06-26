#ifndef SORTEDDATATABLE_H
#define SORTEDDATATABLE_H

#include <set>
#include "include/datasample.h"

class SortedDataTable
{
public:
    SortedDataTable();
    SortedDataTable(bool allowDuplicates);

    void addSample(const DataSample &sample);

    std::multiset<DataSample>::const_iterator cbegin() const;
    std::multiset<DataSample>::const_iterator cend() const;

    unsigned int getXDimension() const;
    unsigned int getYDimension() const;
    unsigned int getNumSamples() const;

    /*
     * Backwards compatibility
     */
    void getBackwardsCompatibleGrid(std::vector< std::vector<double> > &backwardsCompatibleGrid) const; // Creates backwards compatible grid if it doesn't exist
    void getBackwardsCompatibleX(std::vector< std::vector<double> > &backwardsCompatibleX) const;
    void getBackwardsCompatibleY(std::vector< std::vector<double> > &backwardsCompatibleY) const;

    std::vector< std::vector<double> > transposeTable(const std::vector< std::vector<double> > &table) const;
    std::vector< std::vector<double> > getTransposedTableX() const;
    std::vector< std::vector<double> > getTransposedTableY() const;

    /*
     * Debug
     */
    void printSamples() const;  // Print point and values for all samples
    void printGrid() const;     // Print the grid (that we know of so far)

    bool isGridComplete() const;

private:
    bool allowDuplicates;
    unsigned int numDuplicates;
    unsigned int xDim;
    unsigned int yDim;

    std::multiset<DataSample> samples;
    std::vector< std::set<double> > grid;

    void initDataStructures(); // Initialise grid to be a std::vector of xDim std::sets
    unsigned int getNumSamplesRequired() const;

    void recordGridPoint(const DataSample &sample);

    // Used by functions that require the grid to be complete before they start their operation
    // This function prints a message and exits the program if the grid is not complete.
    void gridCompleteGuard() const;
};

#endif // SORTEDDATATABLE_H
