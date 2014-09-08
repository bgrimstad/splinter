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


#ifndef RADIALBASISFUNCTION_H
#define RADIALBASISFUNCTION_H

#include "datatable.h"
#include "include/spline.h"
#include "memory"

namespace MultivariateSplines
{

enum class RadialBasisFunctionType
{
    MULTIQUADRIC,
    INVERSE_QUADRIC,
    INVERSE_MULTIQUADRIC,
    THIN_PLATE_SPLINE,
    GAUSSIAN
};

/*
 * Base class for radial basis functions.
 */
class RadialBasisFunction
{
public:
    RadialBasisFunction() : e(1.0) {}
    RadialBasisFunction(double e) : e(e) {}
    virtual double eval(double r) const = 0;
    // virtual double grad(double r) const = 0; // TODO: implement gradient
protected:
    double e;
};

class ThinPlateSpline : public RadialBasisFunction
{
public:
    double eval(double r) const
    {
        return (r<=0.0) ? 0.0 : r*r*std::log(r);
    }
};

class Multiquadric : public RadialBasisFunction
{
public:
    double eval(double r) const
    {
        return std::sqrt(1.0 + e*e*r*r);
    }
};

class InverseMultiquadric : public RadialBasisFunction
{
public:
    double eval(double r) const
    {
        return 1.0/std::sqrt(1.0 + e*e*r*r);
    }
};

class InverseQuadric : public RadialBasisFunction
{
public:
    double eval(double r) const
    {
        return 1.0/(1.0 + e*e*r*r);
    }
};

class Gaussian : public RadialBasisFunction
{
public:
    double eval(double r) const
    {
        return std::exp(-e*e*r*r);
    }
};

/*
 * Class for radial basis function splines.
 * The RBF splines support scattered sampling, but their construction require
 * the solution of a fairly ill-conditioned linear system. This drawback may be
 * alleviated by applying a pre-conditioner to the linear system.
 */
class RBFSpline : public Spline
{
public:

    RBFSpline(DataTable &samples, RadialBasisFunctionType type);
    RBFSpline(DataTable &samples, RadialBasisFunctionType type, bool normalized);

    virtual RBFSpline* clone() const { return new RBFSpline(*this); }

    double eval(DenseVector &x) const;
    double eval(std::vector<double> &x) const;

    DenseMatrix evalJacobian(DenseVector &x) const {}; // TODO: implement via RBF_fn
    DenseMatrix evalHessian(DenseVector &x) const {}; // TODO: implement via RBF_fn
    //    std::vector<double> getDomainUpperBound() const;
    //    std::vector<double> getDomainLowerBound() const;

    unsigned int getNumVariables() const { return dim; }

private:

    const DataTable samples;
    bool normalized, precondition;
    unsigned int dim, numSamples;

    std::shared_ptr<RadialBasisFunction> fn;

    DenseMatrix weights;

    DenseMatrix computePreconditionMatrix() const;

    double dist(const std::vector<double> x, const std::vector<double> y) const;
    double dist(const DataSample &x, const DataSample &y) const;
    bool dist_sort(const DataSample &x, const DataSample &y) const;

};

/*
 * Helper class for pre-condition matrix calculation
 */
class Point
{
public:
    Point(std::vector<double> p) : p(p),i(0) {}

    bool operator<(const Point &rhs) const
    {
        // Compare dist to origin
        assert(getDim() == rhs.getDim());
        if(dist() < rhs.dist())
            return true;
        return false;
    }

    bool operator==(const Point &rhs) const
    {
        if (getDim() != rhs.getDim())
            return false;

        for (unsigned int i = 0; i < getDim(); i++)
            if (p.at(i) != rhs.getPoint().at(i))
                return false;

        return true;
    }

    Point& operator+=(const Point &rhs)
    {
        assert(getDim() == rhs.getDim());
        for (unsigned int i = 0; i < p.size(); i++)
            p.at(i) += rhs.getPoint().at(i);
        return (*this);
    }

    Point& operator-=(const Point &rhs)
    {
        assert(getDim() == rhs.getDim());
        for (unsigned int i = 0; i < p.size(); i++)
            p.at(i) -= rhs.getPoint().at(i);
        return (*this);
    }

    Point operator+(const Point &rhs) const
    {
        assert(getDim() == rhs.getDim());
        Point temp(p);
        temp += rhs;
        return temp;
    }

    Point operator-(const Point &rhs) const
    {
        assert(getDim() == rhs.getDim());
        Point temp(p);
        temp -= rhs;
        return temp;
    }

    double dist(const Point &rhs) const
    {
        assert(getDim() == rhs.getDim());
        double sum = 0;
        for(unsigned int i=0; i<getDim(); i++)
            sum += (p.at(i)-rhs.getPoint().at(i))*(p.at(i)-rhs.getPoint().at(i));
        return std::sqrt(sum);
    }

    double dist() const
    {
        std::vector<double> origin(getDim(),0.0);
        return dist(Point(origin));
    }

    //friend std::ostream &operator<<(std::ostream &outputStream, const DataSample &sample);

    std::vector<double> getPoint() const { return p; }
    unsigned int getDim() const { return p.size(); }
    unsigned int getIndex() const { return i; }
    void setIndex(unsigned int index) { i=index; }

private:
    std::vector<double> p;
    unsigned int i;
};

} // namespace MultivariateSplines

#endif // RADIALBASISFUNCTION_H
