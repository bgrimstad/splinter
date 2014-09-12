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


#ifndef MS_BSPLINEBASIS1D_H
#define MS_BSPLINEBASIS1D_H

#include "generaldefinitions.h"

namespace MultivariateSplines
{

enum class KnotVectorType
{
    EXPLICIT,   // Knot sequence is explicitly given (it should be regular)
    REGULAR,    // p+1-regular knots (knots are added if necessary)
    FREE,       // Knots for cubic spline interpolation with free end conditions
    EQUIDISTANT // Equidistant p+1-regular knot sequence for all degrees
};

class Basis1D
{
public:
    Basis1D(std::vector<double> &x, unsigned int degree);
    Basis1D(std::vector<double> &x, unsigned int degree, KnotVectorType knotVectorType);

    // Evaluation of basis functions
    SparseVector evaluate(const double x) const;
    SparseVector evaluateDerivative(double x, int r) const;
    DenseVector evaluateFirstDerivative(double x) const; // Depricated

    // Knot vector related
    bool refineKnots(SparseMatrix &A);
    bool insertKnots(SparseMatrix &A, double tau, unsigned int multiplicity = 1);
    // bool insertKnots(SparseMatrix &A, std::vector<tuple<double,int>> newKnots); // Add knots at several locations
    unsigned int knotMultiplicity(const double &tau) const; // Returns the number of repetitions of tau in the knot vector

    // Support related
    bool insideSupport(const double &x) const;
    bool reduceSupport(double lb, double ub, SparseMatrix &A);

    // Getters
    std::vector<double> getKnotVector() const { return knots; }
    unsigned int getBasisDegree() const { return degree; }
    double getKnotValue(unsigned int index) const;
    unsigned int numBasisFunctions() const;
    unsigned int numBasisFunctionsTarget() const;

    // Index getters
    std::vector<int> indexSupportedBasisfunctions(double x) const;
    int indexHalfopenInterval(double x) const;
    unsigned int indexLongestInterval() const;
    unsigned int indexLongestInterval(const std::vector<double> &vec) const;

private:

    // DeBoorCox algorithm for evaluating basis functions
    double deBoorCox(double x, int i, int k) const;
    double deBoorCoxCoeff(double x, double x_min, double x_max) const;

    // Builds basis matrix for alternative evaluation of basis functions
    SparseMatrix buildBasisMatrix(double x, int u, int k, bool diff = false) const;

    // Builds knot insertion matrix
    bool buildKnotInsertionMatrix(SparseMatrix &A, const std::vector<double> &refinedKnots) const;

    // Helper functions
    bool inHalfopenInterval(double x, double x_min, double x_max) const;
    std::vector<double> linspace(double start, double stop, unsigned int points) const;

    // Knot vector related
    bool isKnotVectorRegular() const;
    bool isKnotVectorRegular(const std::vector<double> &vec) const;
    bool isRefinement(const std::vector<double> &refinedKnots) const;
    std::vector<double> knotVectorEquidistant(std::vector<double> &X) const;
    std::vector<double> knotVectorRegular(std::vector<double>& X) const;
    std::vector<double> knotVectorFree(std::vector<double>& X) const;

    // Member variables
    unsigned int degree;
    std::vector<double> knots;
    unsigned int targetNumBasisfunctions;
};

} // namespace MultivariateSplines

#endif // MS_BSPLINEBASIS1D_H
