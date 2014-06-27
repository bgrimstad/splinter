#include "include/datasample.h"

DataSample::DataSample(double x, double y)
{
    setData(std::vector<double>(1, x), std::vector<double>(1, y));
}

DataSample::DataSample(std::vector<double> x, double y)
{
    setData(x, std::vector<double>(1, y));
}

DataSample::DataSample(std::vector<double> x, std::vector<double> y)
{
    setData(x, y);
}

DataSample::DataSample(DenseVector x, double y)
{
    std::vector<double> newX, newY;

    for (int i = 0; i < x.size(); i++)
    {
        newX.push_back(x(i));
    }

    newY.push_back(y);

    setData(newX, newY);
}

DataSample::DataSample(DenseVector x, DenseVector y)
{
    std::vector<double> newX, newY;

    for (int i = 0; i < x.size(); i++)
    {
        newX.push_back(x(i));
    }
    for (int i = 0; i < y.size(); i++)
    {
        newY.push_back(y(i));
    }

    setData(newX, newY);
}

void DataSample::setData(const std::vector<double> &x, const std::vector<double> &y)
{
    this->x = x;
    this->y = y;
}

bool DataSample::operator<(const DataSample &rhs) const
{
    assert(this->getDimX() == rhs.getDimX());

    for(unsigned int i = 0; i < this->getDimX(); i++)
    {
        if(x.at(i) < rhs.getX().at(i))
            return true;
        else if(x.at(i) > rhs.getX().at(i))
            return false;
    }

    return false;
}

std::ostream &operator<<(std::ostream &outputStream, const DataSample &sample)
{
    outputStream << "Sample: (";

    bool firstLoop = true;
    for(auto &coordinate : sample.getX())
    {
        if(!firstLoop)
            outputStream << ", ";

        outputStream << coordinate;
        firstLoop = false;
    }

    outputStream << ") = (";

    firstLoop = true;
    for(auto &value : sample.getY())
    {
        if(!firstLoop)
            outputStream << ", ";

        outputStream << value;
        firstLoop = false;
    }

    outputStream << ")";

    return outputStream;
}
