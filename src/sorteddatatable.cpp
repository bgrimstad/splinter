#include "sorteddatatable.h"

SortedDataTable::SortedDataTable()
    : SortedDataTable(false)
{
}

SortedDataTable::SortedDataTable(bool allowDuplicates)
    : allowDuplicates(allowDuplicates), numDuplicates(0), xDim(0)
{
}

void SortedDataTable::addSample(const DataSample &sample)
{
    if(getNumSamples() == 0)
    {
        xDim = sample.getPoint().size();
        yDim = sample.getValue().size();
        initDataStructures();
    }

    assert(sample.getDimX() == xDim); // All points must have the same dimension
    assert(sample.getDimY() == yDim); // All points must have the same dimension

    // Check if the sample has been added already
    if(samples.count(sample) > 0)
    {
        if(!allowDuplicates)
        {
            std::cout << "Discarding duplicate sample because allowDuplicates is false!" << std::endl;
            std::cout << "Initialise with SortedDataTable(true) to set it to true." << std::endl;
            return;
        }

        numDuplicates++;
    }

    samples.insert(sample);

    recordGridPoint(sample);
}

void SortedDataTable::recordGridPoint(const DataSample &sample)
{
    for(unsigned int i = 0; i < getXDimension(); i++)
    {
        grid.at(i).insert(sample.getPoint().at(i));
    }
}

unsigned int SortedDataTable::getNumSamplesRequired() const
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

bool SortedDataTable::isGridComplete() const
{
    return samples.size() > 0 && samples.size() - numDuplicates == getNumSamplesRequired();
}

void SortedDataTable::initDataStructures()
{
    for(unsigned int i = 0; i < getXDimension(); i++)
    {
        grid.push_back(std::set<double>());
    }
}

unsigned int SortedDataTable::getXDimension() const
{
    return xDim;
}

unsigned int SortedDataTable::getYDimension() const
{
    return yDim;
}

unsigned int SortedDataTable::getNumSamples() const
{
    return samples.size();
}

std::multiset<DataSample>::const_iterator SortedDataTable::cbegin() const
{
    gridCompleteGuard();

    return samples.cbegin();
}

std::multiset<DataSample>::const_iterator SortedDataTable::cend() const
{
    gridCompleteGuard();

    return samples.cend();
}

void SortedDataTable::gridCompleteGuard() const
{
    if(!isGridComplete())
    {
        std::cout << "The grid is not complete yet!" << std::endl;
        exit(1);
    }
}

/***************************
 * Backwards compatibility *
 ***************************/

// = [x1, x2, x3], where x1 = [...] holds unique values of x1 variable
void SortedDataTable::getBackwardsCompatibleGrid(std::vector< std::vector<double> > &backwardsCompatibleGrid) const
{
    gridCompleteGuard();

    backwardsCompatibleGrid.clear();

    for(auto &variable : grid)
    {
        backwardsCompatibleGrid.push_back(std::vector<double>());
        for(auto &value : variable)
        {
            backwardsCompatibleGrid.at(backwardsCompatibleGrid.size() - 1).push_back(value);
        }
    }
}

// = [x1, x2, x3], where x1 = [x11, x12, ...] holds the variables of point x1.
//   1 <= size(x1) = size(x2) = ... = size(xn) < m
void SortedDataTable::getBackwardsCompatibleX(std::vector< std::vector<double> > &backwardsCompatibleX) const
{
    gridCompleteGuard();

    backwardsCompatibleX.clear();

    for(auto &sample : samples)
    {
        backwardsCompatibleX.push_back(sample.getPoint());
    }
}

// = [y1, y2, y3], where y1 = [y11, y12, ...] holds the variables of value y1.
//   1 <= size(yn) < m
void SortedDataTable::getBackwardsCompatibleY(std::vector< std::vector<double> > &backwardsCompatibleY) const
{
    gridCompleteGuard();

    backwardsCompatibleY.clear();

    for(auto &sample : samples)
    {
        backwardsCompatibleY.push_back(sample.getValue());
    }
}

// Copy and transpose table: turn columns to rows and rows to columns
std::vector< std::vector<double> > SortedDataTable::transposeTable(const std::vector< std::vector<double> > &table) const
{
    if (table.size() > 0)
    {
        // Assumes square tables!
        std::vector< std::vector<double> > transp(table.at(0).size());

        for (unsigned int i = 0; i < table.size(); i++)
        {
            for (unsigned int j = 0; j < table.at(i).size(); j++)
            {
                transp.at(j).push_back(table.at(i).at(j));
            }
        }

        return transp;
    }

    return table;
}

std::vector< std::vector<double> > SortedDataTable::getTransposedTableX() const
{
    std::vector< std::vector<double> > backwardsCompatibleX;
    getBackwardsCompatibleX(backwardsCompatibleX);
    return transposeTable(backwardsCompatibleX);
}

std::vector< std::vector<double> > SortedDataTable::getTransposedTableY() const
{
    std::vector< std::vector<double> > backwardsCompatibleY;
    getBackwardsCompatibleY(backwardsCompatibleY);
    return transposeTable(backwardsCompatibleY);
}

/**************
 * Debug code *
 **************/
void SortedDataTable::printSamples() const
{
    for(auto &sample : samples)
    {
        std::cout << sample << std::endl;
    }
}

void SortedDataTable::printGrid() const
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
