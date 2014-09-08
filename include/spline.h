#ifndef SPLINE_H
#define SPLINE_H

#include "generaldefinitions.h"

namespace MultivariateSplines
{

/*
 * Interface for spline classes
 */
class Spline
{
public:
    Spline() {}
    virtual ~Spline() {}

    /*
     * Returns the spline value at x
     */
    virtual double eval(DenseVector &x) const = 0;

    /*
     * Returns the (1 x numVariables) Jacobian evaluated at x
     */
    virtual DenseMatrix evalJacobian(DenseVector &x) const = 0;

    /*
     * Returns the (numVariables x numVariables) Hessian evaluated at x
     */
    virtual DenseMatrix evalHessian(DenseVector &x) const = 0;

};

} // namespace MultivariateSplines

#endif // SPLINE_H
