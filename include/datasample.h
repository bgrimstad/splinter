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


#ifndef MS_DATASAMPLE_H
#define MS_DATASAMPLE_H

#include "generaldefinitions.h"

namespace MultivariateSplines
{

/* Class representing a data sample (x,y)
 * where y is the value obtained by sampling
 * at a point x. Both x and y may be scalar
 * or vectors.
*/
class DataSample
{
public:
    DataSample(double x, double y);
    DataSample(std::vector<double> x, double y);
    DataSample(DenseVector x, double y);

    bool operator<(const DataSample &rhs) const; // Returns false if the two are equal
    friend std::ostream &operator<<(std::ostream &outputStream, const DataSample &sample);

    std::vector<double> getX() const { return x; }
    double getY() const { return y; }
    unsigned int getDimX() const { return x.size(); }

private:
    std::vector<double> x;
    double y;

    void setData(const std::vector<double> &x, const double &y);
};

} // namespace MultivariateSplines

#endif // MS_DATASAMPLE_H
