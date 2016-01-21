/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "rbfnetworkbuilder.h"

namespace SPLINTER
{

RBFNetwork RBFNetwork::Builder::build() const
{
    RBFNetwork network;

    std::shared_ptr<RBF> _fn;
    if (_type == RBFType::THIN_PLATE_SPLINE)
    {
        _fn = std::shared_ptr<RBF>(new ThinPlateSpline());
    }
    else if (_type == RBFType::MULTIQUADRIC)
    {
        _fn = std::shared_ptr<RBF>(new Multiquadric());
    }
    else if (_type == RBFType::INVERSE_QUADRIC)
    {
        _fn = std::shared_ptr<RBF>(new InverseQuadric());
    }
    else if (_type == RBFType::INVERSE_MULTIQUADRIC)
    {
        _fn = std::shared_ptr<RBF>(new InverseMultiquadric());
    }
    else if (_type == RBFType::GAUSSIAN)
    {
        _fn = std::shared_ptr<RBF>(new Gaussian());
    }
    else
    {
        _fn = std::shared_ptr<RBF>(new ThinPlateSpline());
    }

    /* Want to solve the linear system A*w = b,
     * where w is the vector of weights.
     * NOTE: the system is dense and by default badly conditioned.
     * It should be solved by a specialized solver such as GMRES
     * with preconditioning (e.g. ACBF) as in Matlab.
     * NOTE: Consider trying the Łukaszyk–Karmowski metric (for two variables)
     */
    unsigned int numSamples = _data.getNumSamples();
    DenseMatrix A; A.setZero(numSamples, numSamples);
    DenseVector b; b.setZero(numSamples);

    int i=0;
    for (auto it1 = _data.cbegin(); it1 != _data.cend(); ++it1, ++i)
    {
        double sum = 0;
        int j=0;
        for (auto it2 = _data.cbegin(); it2 != _data.cend(); ++it2, ++j)
        {
            double val = _fn->eval(dist(*it1, *it2));
            if (val != 0)
            {
                //A.insert(i,j) = val;
                A(i,j) = val;
                sum += val;
            }
        }

        double y = it1->getY();
        if (_normalized) b(i) = sum*y;
        else b(i) = y;
    }

    if (_precondition)
    {
        // Calculate precondition matrix P
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

    // Solve for coefficients/weights
    auto coefficients = svd.solve(b);

#ifndef NDEBUG
    // Compute error. If it is used later on, move this statement above the NDEBUG
    double err = (A * coefficients - b).norm() / b.norm();
    std::cout << "Error: " << std::setprecision(10) << err << std::endl;
#endif // NDEBUG

//    // Alternative solver
//    DenseQR s;
//    bool success = s.solve(A,b,coefficients);
//    assert(success);

    // NOTE: Tried using experimental GMRES solver in Eigen, but it did not work very well.

    network.numVariables = _data.getNumVariables();
    network.coefficients = coefficients;
    network.samples = _data;
    network.normalized = _normalized;
    network.type = _type;
    network.fn = _fn;

    return network;
}


/*
 * Calculate precondition matrix
 * TODO: implement
 */
DenseMatrix RBFNetwork::Builder::computePreconditionMatrix() const
{
    return DenseMatrix::Zero(_data.getNumSamples(), _data.getNumSamples());
}

} // namespace SPLINTER