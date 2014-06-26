#ifndef DATASAMPLE_H
#define DATASAMPLE_H

#include "generaldefinitions.h"

class DataSample
{
public:
    DataSample(double point, double value);
    DataSample(std::vector<double> point, double value);
    DataSample(std::vector<double> point, std::vector<double> value);
    DataSample(DenseVector point, DenseVector value);

    bool operator<(const DataSample &rhs) const; // Returns false if the two are equal
    friend std::ostream &operator<<(std::ostream &outputStream, const DataSample &sample);

    std::vector<double> getPoint() const;
    std::vector<double> getValue() const;
    unsigned int getDimX() const;
    unsigned int getDimY() const;

private:
    std::vector<double> point;
    std::vector<double> value;

    void setData(const std::vector<double> &point, const std::vector<double> &value);
};

#endif // DATASAMPLE_H
