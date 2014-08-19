#ifndef SORTEDDATATABLE_H
#define SORTEDDATATABLE_H

#include <set>
#include "include/datasample.h"
#include <string>

#define SAVE_DOUBLE_PRECISION 9

class SortedDataTable
{
public:
    SortedDataTable();
    SortedDataTable(bool allowDuplicates);

    void addSample(const DataSample &sample);
    void addSample(double x, double y);
    void addSample(std::vector<double> x, double y);
    void addSample(std::vector<double> x, std::vector<double> y);
    void addSample(DenseVector x, double y);
    void addSample(DenseVector x, DenseVector y);

    std::multiset<DataSample>::const_iterator cbegin() const;
    std::multiset<DataSample>::const_iterator cend() const;

    unsigned int getDimX() const;
    unsigned int getDimY() const;
    unsigned int getNumSamples() const;

    void saveDataTable(std::string fileName) const; // Throws std::ios_base::failure
    void loadDataTable(std::string fileName);       // Throws std::ios_base::failure

    /*
     * Backwards compatibility
     */
    std::vector< std::vector<double> > getBackwardsCompatibleGrid() const; // Creates backwards compatible grid if it doesn't exist
    std::vector< std::vector<double> > getBackwardsCompatibleX() const;
    std::vector< std::vector<double> > getBackwardsCompatibleY() const;

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
