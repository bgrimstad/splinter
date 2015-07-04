/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_RADIALBASISFUNCTION_H
#define SPLINTER_RADIALBASISFUNCTION_H

#include "datatable.h"
#include "approximant.h"
#include "radialbasisfunctionterm.h"
#include "memory"

namespace SPLINTER
{

/*
 * Class for radial basis function splines.
 * The RBF splines support scattered sampling, but their construction require
 * the solution of a fairly ill-conditioned linear system. This drawback may be
 * alleviated by applying a pre-conditioner to the linear system.
 */
class API RadialBasisFunction : public Approximant
{
public:
    RadialBasisFunction(const char *filename);
    RadialBasisFunction(const std::string filename);
    RadialBasisFunction(const DataTable &samples, RadialBasisFunctionType type);
    RadialBasisFunction(const DataTable &samples, RadialBasisFunctionType type, bool normalized);

    virtual RadialBasisFunction* clone() const { return new RadialBasisFunction(*this); }

    double eval(DenseVector x) const;
    double eval(std::vector<double> x) const;

    DenseMatrix evalJacobian(DenseVector x) const { return DenseMatrix(); }; // TODO: implement
    DenseMatrix evalHessian(DenseVector x) const { return DenseMatrix(); }; // TODO: implement
    //    std::vector<double> getDomainUpperBound() const;
    //    std::vector<double> getDomainLowerBound() const;

    unsigned int getNumVariables() const override { return dim; }

    void save(const std::string fileName) const override { throw Exception("RadialBasisFunction::save: not implemented."); }
    void load(const std::string fileName) override { throw Exception("RadialBasisFunction::load: not implemented."); }

private:

    const DataTable samples;
    bool normalized, precondition;
    unsigned int dim, numSamples;

    std::shared_ptr<RadialBasisFunctionTerm> fn;

    DenseMatrix weights;

    DenseMatrix computePreconditionMatrix() const;

    double dist(std::vector<double> x, std::vector<double> y) const;
    double dist(DataSample x, DataSample y) const;
    bool dist_sort(DataSample x, DataSample y) const;

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
        if (dist() < rhs.dist())
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
        for (unsigned int i=0; i<getDim(); i++)
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

} // namespace SPLINTER

#endif // SPLINTER_RADIALBASISFUNCTION_H
