/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINEBASIS1D_H
#define SPLINTER_BSPLINEBASIS1D_H

#include "generaldefinitions.h"

namespace SPLINTER
{

class BSplineBasis1D
{
public:
    BSplineBasis1D();
    BSplineBasis1D(std::vector<double> &x, unsigned int degree);
    BSplineBasis1D(std::vector<double> &x, unsigned int degree, bool explicitKnots);
    BSplineBasis1D(unsigned int degree, std::vector<double> &knots, unsigned int targetNumBasisFunctions);

    // Evaluation of basis functions
    SparseVector evaluate(double x) const;
    SparseVector evaluateDerivative(double x, int r) const;
    SparseVector evaluateFirstDerivative(double x) const; // Depricated

    // Knot vector related
    SparseMatrix refineKnots();
    SparseMatrix refineKnotsLocally(double x);
    SparseMatrix decomposeToBezierForm();
    SparseMatrix insertKnots(double tau, unsigned int multiplicity = 1);
    // bool insertKnots(SparseMatrix &A, std::vector<tuple<double,int>> newKnots); // Add knots at several locations
    unsigned int knotMultiplicity(double tau) const; // Returns the number of repetitions of tau in the knot vector

    // Support related
    void supportHack(double &x) const;
    bool insideSupport(double x) const;
    SparseMatrix reduceSupport(double lb, double ub);

    // Getters
    std::vector<double> getKnotVector() const { return knots; }
    unsigned int getBasisDegree() const { return degree; }
    double getKnotValue(unsigned int index) const;
    unsigned int getNumBasisFunctions() const;
    unsigned int getNumBasisFunctionsTarget() const;
    unsigned int getDegree() const;
    std::vector<double> getKnots() const;
    unsigned int getTargetNumBasisFunctions() const;

    // Index getters
    std::vector<int> indexSupportedBasisfunctions(double x) const;
    int indexHalfopenInterval(double x) const;
    unsigned int indexLongestInterval() const;
    unsigned int indexLongestInterval(const std::vector<double> &vec) const;

    // Setters
    void setNumBasisFunctionsTarget(unsigned int target)
    {
        targetNumBasisfunctions = std::max(degree+1, target);
    }

private:

    // DeBoorCox algorithm for evaluating basis functions
    double deBoorCox(double x, int i, int k) const;
    double deBoorCoxCoeff(double x, double x_min, double x_max) const;

    // Builds basis matrix for alternative evaluation of basis functions
    SparseMatrix buildBasisMatrix(double x, unsigned int u, unsigned int k, bool diff = false) const;

    // Builds knot insertion matrix
    SparseMatrix buildKnotInsertionMatrix(const std::vector<double> &refinedKnots) const;

    // Helper functions
    bool inHalfopenInterval(double x, double x_min, double x_max) const;

    // Knot vector related
    bool isKnotVectorRegular() const;
    bool isKnotVectorRegular(const std::vector<double> &vec) const;
    bool isRefinement(const std::vector<double> &refinedKnots) const;
    std::vector<double> knotVectorMovingAverage(std::vector<double> &vec) const;

    // Member variables
    unsigned int degree;
    std::vector<double> knots;
    unsigned int targetNumBasisfunctions;

    friend class Serializer;
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINEBASIS1D_H
