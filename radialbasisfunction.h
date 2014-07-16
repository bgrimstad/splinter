#ifndef RADIALBASISFUNCTION_H
#define RADIALBASISFUNCTION_H

#include "memory"
#include "sorteddatatable.h"

enum class RadialBasisFunctionType
{
    MULTIQUADRIC,
    INVERSE_QUADRIC,
    INVERSE_MULTIQUADRIC,
    THIN_PLATE_SPLINE,
    GAUSSIAN
};

class RBF_fn
{
public:
    RBF_fn() : e(1.0) {}
    RBF_fn(double e) : e(e) {}
    virtual double eval(double r) const = 0;
protected:
    double e;
};

class ThinPlateSpline : public RBF_fn
{
public:
    double eval(double r) const
    {
        return (r==0.0) ? 0.0 : r*r*std::log(r);
    }
};

class Multiquadric : public RBF_fn
{
public:
    double eval(double r) const
    {
        return std::sqrt(1 + e*e*r*r);
    }
};

class InverseQuadric : public RBF_fn
{
public:
    double eval(double r) const
    {
        return 1.0/(1.0 + e*e*r*r);
    }
};

class Gaussian : public RBF_fn
{
public:
    double eval(double r) const
    {
        return std::exp(-e*e*r*r);
    }
};

class InverseMultiquadric : public RBF_fn
{
public:
    double eval(double r) const
    {
        return 1.0/std::sqrt(1 + e*e*r*r);
    }
};

class RadialBasisFunction
{
public:

    RadialBasisFunction(SortedDataTable &samples, RadialBasisFunctionType type);
    RadialBasisFunction(SortedDataTable &samples, RadialBasisFunctionType type, bool normalized);

    virtual RadialBasisFunction* clone() const { return new RadialBasisFunction(*this); }

    double evaluate(DenseVector &x) const;
    double evaluate(std::vector<double> &x) const;
    //    std::vector<double> getDomainUpperBound() const;
    //    std::vector<double> getDomainLowerBound() const;

    int dimX() const;

private:

    const SortedDataTable samples;
    bool normalized;
    unsigned int dim, numSamples;

    std::shared_ptr<RBF_fn> fn;

    DenseMatrix weights;

    double dist(std::vector<double> x, std::vector<double> y) const;
};

#endif // RADIALBASISFUNCTION_H
