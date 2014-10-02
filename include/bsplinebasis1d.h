/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
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

class BSplineBasis1D
{
public:
    BSplineBasis1D(std::vector<double> &x, unsigned int degree);
    BSplineBasis1D(std::vector<double> &x, unsigned int degree, KnotVectorType knotVectorType);

    // Evaluation of basis functions
    SparseVector evaluate(double x) const;
    SparseVector evaluateDerivative(double x, int r) const;
    DenseVector evaluateFirstDerivative(double x) const; // Depricated

    // Knot vector related
    bool refineKnots(SparseMatrix &A);
    bool insertKnots(SparseMatrix &A, double tau, unsigned int multiplicity = 1);
    // bool insertKnots(SparseMatrix &A, std::vector<tuple<double,int>> newKnots); // Add knots at several locations
    unsigned int knotMultiplicity(double tau) const; // Returns the number of repetitions of tau in the knot vector

    // Support related
    void supportHack(double &x) const;
    bool insideSupport(double x) const;
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
    SparseMatrix buildBasisMatrix(double x, unsigned int u, unsigned int k, bool diff = false) const;

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
