#ifndef PSPLINE_H
#define PSPLINE_H

#include "bspline.h"

namespace MultivariateSplines
{

/*
 * The P-Spline is a smooting spline which relaxes the interpolation constraints on the control points to allow smoother spline curves.
 * It minimizes objective which penalizes both deviation (for interpolation) and second derivative (for smoothing).
 * It inherits all properties of the B-spline - the only difference lies in the calculation of the control points.
 */
class PSpline : public BSpline
{
public:
    PSpline(DataTable &samples);
    PSpline(DataTable &samples, double lambda);

private:

    // Smoothing parameter (usually set to a small number; default 0.03)
    double lambda;

    // P-spline control point calculation
    void computeControlPoints(const DataTable &samples, std::vector< std::vector<double> > &X, std::vector< std::vector<double> > &Y);
    void controlPointEquationLHS(const DataTable &samples, std::vector< std::vector<double> > &X, SparseMatrix &L, SparseMatrix &B, SparseMatrix &W, double lambda);
    void controlPointEquationRHS(std::vector< std::vector<double> > &Y, DenseMatrix &R, SparseMatrix &B, SparseMatrix &W);
    void getSecondOrderFiniteDifferenceMatrix(SparseMatrix &D);
    void getBasisFunctionMatrix(std::vector< std::vector<double> > &X, SparseMatrix &B);

};

} // namespace MultivariateSplines

#endif // PSPLINE_H
