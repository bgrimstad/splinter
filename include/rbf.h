/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_RBF_H
#define SPLINTER_RBF_H

#include "datatable.h"
#include "definitions.h"
#include "memory"

namespace SPLINTER
{

enum class RBFType
{
    MULTIQUADRIC,
    INVERSE_QUADRIC,
    INVERSE_MULTIQUADRIC,
    THIN_PLATE_SPLINE,
    GAUSSIAN
};

/*
 * Base class for radial basis functions.
 * TODO: this class should inherit from Function
 */
class RBF
{
public:
    RBF() : e(1.0) {}
    RBF(double e) : e(e) {}
    virtual double eval(double r) const = 0;
    virtual double evalDerivative(double r) const = 0;
    virtual ~RBF() {}

protected:
    double e;
};

class ThinPlateSpline : public RBF
{
public:
    double eval(double r) const
    {
        return (r<=0.0) ? 0.0 : r*r*std::log(r);
    }
    double evalDerivative(double r) const
    {
        return (r<=0.0) ? 0.0 : r*(2*std::log(r) + 1);
    }
};

class Multiquadric : public RBF
{
public:
    double eval(double r) const
    {
        return std::sqrt(1.0 + e*e*r*r);
    }
    double evalDerivative(double r) const
    {
        return e*e*r/std::sqrt(1 + e*e*r*r);
    }
};

class InverseMultiquadric : public RBF
{
public:
    double eval(double r) const
    {
        return 1.0/std::sqrt(1.0 + e*e*r*r);
    }
    double evalDerivative(double r) const
    {
        return -e*e*r/(std::sqrt(1 + e*e*r*r)*(1 + e*e*r*r));
    }
};

class InverseQuadric : public RBF
{
public:
    double eval(double r) const
    {
        return 1.0/(1.0 + e*e*r*r);
    }
    double evalDerivative(double r) const
    {
        return -2*e*e*r/((1 + e*e*r*r)*(1 + e*e*r*r));
    }
};

class Gaussian : public RBF
{
public:
    double eval(double r) const
    {
        return std::exp(-e*e*r*r);
    }
    double evalDerivative(double r) const
    {
        return -2*e*e*r*std::exp(-e*e*r*r);
    }
};

} // namespace SPLINTER

#endif // SPLINTER_RBF_H
