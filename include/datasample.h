/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
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

    void setData(const std::vector<double> &x, double y);
};

} // namespace MultivariateSplines

#endif // MS_DATASAMPLE_H
