/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <serializer.h>
#include "rbfapproximant.h"
#include "linearsolvers.h"
#include "Eigen/SVD"

namespace SPLINTER
{

RBFApproximant::RBFApproximant()
    : Function(1)
{
}

RBFApproximant::RBFApproximant(const char *fileName)
    : RBFApproximant(std::string(fileName))
{
}

RBFApproximant::RBFApproximant(const std::string fileName)
    : Function(1)
{
    load(fileName);
}

RBFApproximant::RBFApproximant(const DataTable &samples, RBFType type)
    : RBFApproximant(samples, type, false)
{
}

RBFApproximant::RBFApproximant(const DataTable &samples, RBFType type, bool normalized)
    : Function(samples.getNumVariables()),
      samples(samples),
      type(type),
      normalized(normalized),
      precondition(false),
      numSamples(samples.getNumSamples())
{
    if (type == RBFType::THIN_PLATE_SPLINE)
    {
        fn = std::shared_ptr<RBFTerm>(new ThinPlateSpline());
    }
    else if (type == RBFType::MULTIQUADRIC)
    {
        fn = std::shared_ptr<RBFTerm>(new Multiquadric());
    }
    else if (type == RBFType::INVERSE_QUADRIC)
    {
        fn = std::shared_ptr<RBFTerm>(new InverseQuadric());
    }
    else if (type == RBFType::INVERSE_MULTIQUADRIC)
    {
        fn = std::shared_ptr<RBFTerm>(new InverseMultiquadric());
    }
    else if (type == RBFType::GAUSSIAN)
    {
        fn = std::shared_ptr<RBFTerm>(new Gaussian());
    }
    else
    {
        fn = std::shared_ptr<RBFTerm>(new ThinPlateSpline());
    }

    /* Want to solve the linear system A*w = b,
     * where w is the vector of weights.
     * NOTE: the system is dense and by default badly conditioned.
     * It should be solved by a specialized solver such as GMRES
     * with preconditioning (e.g. ACBF) as in Matlab.
     * NOTE: Consider trying the Łukaszyk–Karmowski metric (for two variables)
     */
    DenseMatrix A; A.setZero(numSamples, numSamples);
    DenseMatrix b; b.setZero(numSamples,1);

    int i=0;
    for (auto it1 = samples.cbegin(); it1 != samples.cend(); ++it1, ++i)
    {
        double sum = 0;
        int j=0;
        for (auto it2 = samples.cbegin(); it2 != samples.cend(); ++it2, ++j)
        {
            double val = fn->eval(dist(*it1, *it2));
            if (val != 0)
            {
                //A.insert(i,j) = val;
                A(i,j) = val;
                sum += val;
            }
        }

        double y = it1->getY();
        if (normalized) b(i) = sum*y;
        else b(i) = y;
    }

    if (precondition)
    {
        // Calcualte precondition matrix P
        DenseMatrix P = computePreconditionMatrix();

        // Preconditioned A and b
        DenseMatrix Ap = P*A;
        DenseMatrix bp = P*b;

        A = Ap;
        b = bp;
    }

    #ifndef NDEBUG
    std::cout << "Computing RBF weights using dense solver." << std::endl;
    #endif // NDEBUG

    // SVD analysis
    Eigen::JacobiSVD<DenseMatrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    auto svals = svd.singularValues();
    double svalmax = svals(0);
    double svalmin = svals(svals.rows()-1);
    double rcondnum = (svalmax <= 0.0 || svalmin <= 0.0) ? 0.0 : svalmin/svalmax;

    #ifndef NDEBUG
    std::cout << "The reciprocal of the condition number is: " << rcondnum << std::endl;
    std::cout << "Largest/smallest singular value: " << svalmax << " / " << svalmin << std::endl;
    #endif // NDEBUG

    // Solve for weights
    weights = svd.solve(b);

    #ifndef NDEBUG
    // Compute error. If it is used later on, move this statement above the NDEBUG
    double err = (A*weights - b).norm() / b.norm();
    std::cout << "Error: " << std::setprecision(10) << err << std::endl;
    #endif // NDEBUG

//    // Alternative solver
//    DenseQR s;
//    bool success = s.solve(A,b,weights);
//    assert(success);

    // NOTE: Tried using experimental GMRES solver in Eigen, but it did not work very well.
}

double RBFApproximant::eval(DenseVector x) const
{
    std::vector<double> y;
    for (int i=0; i<x.rows(); i++)
        y.push_back(x(i));
    return eval(y);
}

double RBFApproximant::eval(std::vector<double> x) const
{
    if (x.size() != numVariables)
        throw Exception("RBFApproximant::eval: Wrong dimension on evaluation point x.");

    double fval, sum = 0, sumw = 0;
    int i = 0;
    for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        fval = fn->eval(dist(x,it->getX()));
        sumw += weights(i)*fval;
        sum += fval;
    }
    return normalized ? sumw/sum : sumw;
}

DenseMatrix RBFApproximant::evalJacobian(DenseVector x) const
{
    std::vector<double> x_vec;
    for (unsigned int i = 0; i<x.size(); i++)
        x_vec.push_back(x(i));

    DenseMatrix jac;
    jac.setZero(1,numVariables);

    for (unsigned int i = 0; i < numVariables; i++)
    {
        double sumw = 0;
        double sumw_d = 0;
        double sum = 0;
        double sum_d = 0;

        int j = 0;
        for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++j)
        {
            // Sample
            auto s_vec = it->getX();

            // Distance from sample
            double r = dist(x_vec, s_vec);
            double ri = x_vec.at(i) - s_vec.at(i);

            // Evaluate RBF and its derivative at r
            double f = fn->eval(r);
            double dfdr = fn->evalDerivative(r);

            sum += f;
            sumw += weights(j)*f;

            // TODO: check if this assumption is correct
            if (r != 0)
            {
                sum_d += dfdr*ri/r;
                sumw_d += weights(j)*dfdr*ri/r;
            }
        }

        if (normalized)
            jac(i) = (sum*sumw_d - sum_d*sumw)/(sum*sum);
        else
            jac(i) = sumw_d;
    }
    return jac;
}

/*
 * Calculate precondition matrix
 * TODO: implement
 */
DenseMatrix RBFApproximant::computePreconditionMatrix() const
{
    return DenseMatrix::Zero(numSamples, numSamples);
}

/*
 * Computes Euclidean distance ||x-y||
 */
double RBFApproximant::dist(std::vector<double> x, std::vector<double> y) const
{
    if (x.size() != y.size())
        throw Exception("RBFApproximant::dist: Cannot measure distance between two points of different dimension");
    double sum = 0.0;
    for (unsigned int i=0; i<x.size(); i++)
        sum += (x.at(i)-y.at(i))*(x.at(i)-y.at(i));
    return std::sqrt(sum);
}

/*
 * Computes Euclidean distance ||x-y||
 */
double RBFApproximant::dist(DataSample x, DataSample y) const
{
    return dist(x.getX(), y.getX());
}

bool RBFApproximant::dist_sort(DataSample x, DataSample y) const
{
    std::vector<double> zeros(x.getDimX(), 0);
    DataSample origin(zeros, 0.0);
    double x_dist = dist(x, origin);
    double y_dist = dist(y, origin);
    return (x_dist<y_dist);
}

void RBFApproximant::save(const std::string fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void RBFApproximant::load(const std::string fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

const std::string RBFApproximant::getDescription() const {
    std::string description("RadialBasisFunction of type ");
    switch(type) {
    case RBFType::GAUSSIAN:
        description.append("Gaussian");
        break;
    case RBFType::INVERSE_MULTIQUADRIC:
        description.append("Inverse multiquadric");
        break;
    case RBFType::INVERSE_QUADRIC:
        description.append("Inverse quadric");
        break;
    case RBFType::MULTIQUADRIC:
        description.append("Multiquadric");
        break;
    case RBFType::THIN_PLATE_SPLINE:
        description.append("Thin plate spline");
        break;
    default:
        description.append("Error: Unknown!");
        break;
    }

    return description;
}

} // namespace SPLINTER
