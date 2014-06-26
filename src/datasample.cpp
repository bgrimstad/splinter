#include "datasample.h"

DataSample::DataSample(double point, double value)
{
    setData(std::vector<double>(1, point), std::vector<double>(1, value));
}

DataSample::DataSample(std::vector<double> point, double value)
{
    setData(point, std::vector<double>(1, value));
}

DataSample::DataSample(std::vector<double> point, std::vector<double> value)
{
    setData(point, value);
}

DataSample::DataSample(DenseVector point, DenseVector value)
{
    std::vector<double> newPoint, newValue;

    for (int i = 0; i < point.size(); i++)
    {
        newPoint.push_back(point(i));
    }
    for (int i = 0; i < value.size(); i++)
    {
        newValue.push_back(value(i));
    }

    setData(newPoint, newValue);
}

void DataSample::setData(const std::vector<double> &point, const std::vector<double> &value)
{
    this->point = point;
    this->value = value;
}

std::vector<double> DataSample::getPoint() const
{
    return point;
}

std::vector<double> DataSample::getValue() const
{
    return value;
}

unsigned int DataSample::getDimX() const
{
    return point.size();
}

unsigned int DataSample::getDimY() const
{
    return value.size();
}

bool DataSample::operator<(const DataSample &rhs) const
{
    assert(this->getDimX() == rhs.getDimX());

    for(unsigned int i = 0; i < this->getDimX(); i++)
    {
        if(point.at(i) < rhs.getPoint().at(i))
            return true;
        else if(point.at(i) > rhs.getPoint().at(i))
            return false;
    }

    return false;
}

std::ostream &operator<<(std::ostream &outputStream, const DataSample &sample)
{
    outputStream << "Sample: (";

    bool firstLoop = true;
    for(auto &coordinate : sample.getPoint())
    {
        if(!firstLoop)
            outputStream << ", ";

        outputStream << coordinate;
        firstLoop = false;
    }

    outputStream << ") = (";

    firstLoop = true;
    for(auto &value : sample.getValue())
    {
        if(!firstLoop)
            outputStream << ", ";

        outputStream << value;
        firstLoop = false;
    }

    outputStream << ")";

    return outputStream;
}
