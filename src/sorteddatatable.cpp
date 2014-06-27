#include "sorteddatatable.h"

SortedDataTable::SortedDataTable()
    : SortedDataTable(false)
{
}

SortedDataTable::SortedDataTable(bool allowDuplicates)
    : allowDuplicates(allowDuplicates), numDuplicates(0), xDim(0)
{
}

void SortedDataTable::addSample(double x, double y)
{
    addSample(DataSample(x, y));
}

void SortedDataTable::addSample(std::vector<double> x, double y)
{
    addSample(DataSample(x, y));
}

void SortedDataTable::addSample(std::vector<double> x, std::vector<double> y)
{
    addSample(DataSample(x, y));
}

void SortedDataTable::addSample(DenseVector x, double y)
{
    addSample(DataSample(x, y));
}

void SortedDataTable::addSample(DenseVector x, DenseVector y)
{
    addSample(DataSample(x, y));
}

void SortedDataTable::addSample(const DataSample &sample)
{
    if(getNumSamples() == 0)
    {
        xDim = sample.getX().size();
        yDim = sample.getY().size();
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
    for(unsigned int i = 0; i < getDimX(); i++)
    {
        grid.at(i).insert(sample.getX().at(i));
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
    for(unsigned int i = 0; i < getDimX(); i++)
    {
        grid.push_back(std::set<double>());
    }
}

unsigned int SortedDataTable::getDimX() const
{
    return xDim;
}

unsigned int SortedDataTable::getDimY() const
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
std::vector< std::vector<double> > SortedDataTable::getBackwardsCompatibleGrid() const
{
    gridCompleteGuard();

    std::vector< std::vector<double> > backwardsCompatibleGrid;

    for(auto &variable : grid)
    {
        backwardsCompatibleGrid.push_back(std::vector<double>());
        for(auto &value : variable)
        {
            backwardsCompatibleGrid.at(backwardsCompatibleGrid.size() - 1).push_back(value);
        }
    }

    return backwardsCompatibleGrid;
}

// = [x1, x2, x3], where x1 = [x11, x12, ...] holds the variables of point x1.
//   1 <= size(x1) = size(x2) = ... = size(xn) < m
std::vector< std::vector<double> > SortedDataTable::getBackwardsCompatibleX() const
{
    gridCompleteGuard();

    std::vector< std::vector<double> > backwardsCompatibleX;

    for(auto &sample : samples)
    {
        backwardsCompatibleX.push_back(sample.getX());
    }

    return backwardsCompatibleX;
}

// = [y1, y2, y3], where y1 = [y11, y12, ...] holds the variables of value y1.
//   1 <= size(yn) < m
std::vector< std::vector<double> > SortedDataTable::getBackwardsCompatibleY() const
{
    gridCompleteGuard();

    std::vector< std::vector<double> > backwardsCompatibleY;

    for(auto &sample : samples)
    {
        backwardsCompatibleY.push_back(sample.getY());
    }

    return backwardsCompatibleY;
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
    return transposeTable(getBackwardsCompatibleX());
}

std::vector< std::vector<double> > SortedDataTable::getTransposedTableY() const
{
    return transposeTable(getBackwardsCompatibleY());
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
