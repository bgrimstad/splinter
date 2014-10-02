/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "rbfspline.h"
#include "include/linearsolvers.h"
#include "Eigen/SVD"

namespace MultivariateSplines
{

RBFSpline::RBFSpline(const DataTable &samples, RadialBasisFunctionType type)
    : RBFSpline(samples, type, false)
{
}

RBFSpline::RBFSpline(const DataTable &samples, RadialBasisFunctionType type, bool normalized)
    : samples(samples),
      normalized(normalized),
      precondition(false),
      dim(samples.getNumVariables()),
      numSamples(samples.getNumSamples())
{
    if (type == RadialBasisFunctionType::THIN_PLATE_SPLINE)
    {
        fn = std::shared_ptr<RadialBasisFunction>(new ThinPlateSpline());
    }
    else if (type == RadialBasisFunctionType::MULTIQUADRIC)
    {
        fn = std::shared_ptr<RadialBasisFunction>(new Multiquadric());
    }
    else if (type == RadialBasisFunctionType::INVERSE_QUADRIC)
    {
        fn = std::shared_ptr<RadialBasisFunction>(new InverseQuadric());
    }
    else if (type == RadialBasisFunctionType::INVERSE_MULTIQUADRIC)
    {
        fn = std::shared_ptr<RadialBasisFunction>(new InverseMultiquadric());
    }
    else if (type == RadialBasisFunctionType::GAUSSIAN)
    {
        fn = std::shared_ptr<RadialBasisFunction>(new Gaussian());
    }
    else
    {
        fn = std::shared_ptr<RadialBasisFunction>(new ThinPlateSpline());
    }

    /* Want to solve the linear system A*w = b,
     * where w is the vector of weights.
     * NOTE: the system is dense and by default badly conditioned.
     * It should be solved by a specialized solver such as GMRES
     * with preconditioning (e.g. ACBF) as in matlab.
     * NOTE: Consider trying the Łukaszyk–Karmowski metric (for two variables)
     */
    //SparseMatrix A(numSamples,numSamples);
    //A.reserve(numSamples*numSamples);
    DenseMatrix A; A.setZero(numSamples, numSamples);
    DenseMatrix b; b.setZero(numSamples,1);

    int i=0;
    for(auto it1 = samples.cbegin(); it1 != samples.cend(); ++it1, ++i)
    {
        double sum = 0;
        int j=0;
        for(auto it2 = samples.cbegin(); it2 != samples.cend(); ++it2, ++j)
        {
            double val = fn->eval(dist(*it1, *it2));
            if(val != 0)
            {
                //A.insert(i,j) = val;
                A(i,j) = val;
                sum += val;
            }
        }

        double y = it1->getY();
        if(normalized) b(i) = sum*y;
        else b(i) = y;
    }

    //A.makeCompressed();

    if(precondition)
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

double RBFSpline::eval(DenseVector x) const
{
    std::vector<double> y;
    for(int i=0; i<x.rows(); i++)
        y.push_back(x(i));
    return eval(y);
}

double RBFSpline::eval(std::vector<double> x) const
{
    assert(x.size() == dim);
    double fval, sum = 0, sumw = 0;
    int i = 0;
    for(auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        fval = fn->eval(dist(x,it->getX()));
        sumw += weights(i)*fval;
        sum += fval;
    }
    return normalized ? sumw/sum : sumw;
}

/*
 * TODO: test for errors
 */
//DenseMatrix RBFSpline::evalJacobian(DenseVector x) const
//{
//    std::vector<double> x_vec;
//    for(unsigned int i = 0; i<x.size(); i++)
//        x_vec.push_back(x(i));

//    DenseMatrix jac;
//    jac.setZero(1,dim);

//    for(unsigned int i = 0; i < dim; i++)
//    {
//        double sumw = 0;
//        double sumw_d = 0;
//        double sum = 0;
//        double sum_d = 0;

//        int j = 0;
//        for(auto it = samples.cbegin(); it != samples.cend(); ++it, ++j)
//        {
//            // Sample
//            auto s_vec = it->getX();

//            // Distance from sample
//            double r = dist(x_vec, s_vec);
//            double ri = x_vec.at(i) - s_vec.at(i);

//            // Evaluate RBF and its derivative at r
//            double f = fn->eval(r);
//            double dfdr = fn->evalDerivative(r);

//            sum += f;
//            sumw += weights(j)*f;

//            // TODO: check if this assumption is correct
//            if(r != 0)
//            {
//                sum_d += dfdr*ri/r;
//                sumw_d += weights(j)*dfdr*ri/r;
//            }
//        }

//        if(normalized)
//            jac(i) = (sum*sumw_d - sum_d*sumw)/(sum*sum);
//        else
//            jac(i) = sumw_d;
//    }
//    return jac;
//}

/*
 * Calculate precondition matrix
 */
DenseMatrix RBFSpline::computePreconditionMatrix() const
{
    DenseMatrix P;
    P.setZero(numSamples,numSamples);

    // Calculate precondition matrix P based on
    // purely local approximate cardinal basis functions (ACBF)
    int sigma = std::max(1.0, std::floor(0.1*numSamples)); // Local points to consider

    int i=0;
    for(auto it1 = samples.cbegin(); it1 != samples.cend(); ++it1, ++i)
    {
        Point p1(it1->getX());

        // Shift data using p1 as origin
        std::vector<Point> shifted_points;
        int j=0;
        for(auto it2 = samples.cbegin(); it2 != samples.cend(); ++it2, ++j)
        {
            Point p2(it2->getX());
            Point p3(p2-p1);
            p3.setIndex(j); // store index with point
            shifted_points.push_back(p3);
        }
        std::sort(shifted_points.begin(), shifted_points.end());

        // Find sigma closest points to p1
        std::vector<Point> points;
        std::vector<int> indices;
        for(int j=0; j<sigma; j++)
        {
            Point p(shifted_points.at(j));
            indices.push_back(p.getIndex());
            Point p2 = p+p1;
            p2.setIndex(p.getIndex());
            points.push_back(p2); // The resulting point has a different index than that of p
//                cout << p.getIndex() << "/" << p1.getIndex() << "/" << p2.getIndex() << endl;
//                assert(p.getIndex() == p2.getIndex());
        }

        // Add some points far away from p1 and preferably scattered around the domain boundary
        for(int k=0; k<1; k++)
        {
            Point p(shifted_points.at(shifted_points.size()-1-k));
            indices.push_back(p.getIndex());
            Point p2 = p+p1;
            p2.setIndex(p.getIndex());
            points.push_back(p2);
        }

        // Build and solve linear system
        int m = points.size();

        DenseMatrix e; e.setZero(m,1); e(0,0) = 1;
        DenseMatrix B; B.setZero(m,m);
        DenseMatrix w; w.setZero(m,1);

        assert(points.front().getIndex() == i);

        for(int k1=0; k1<m; k1++)
        {
            for(int k2=0; k2<m; k2++)
            {
                Point p = points.at(k1) - points.at(k2);
                B(k1,k2) = fn->eval(p.dist());
            }
        }
//            DenseQR s;
//            bool success = s.solve(B,e,w);
//            assert(success);

        Eigen::JacobiSVD<DenseMatrix> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
        w = svd.solve(e);
        //assert(svd.info() == Eigen::Success);

        for(unsigned int j=0; j<numSamples; j++)
        {
            auto it = find(indices.begin(),indices.end(),j);
            if(it!=indices.end())
            {
                int k = it - indices.begin();

                P(i,j) = w(k,0);
                //cout << "j/k/g " << j << "/" << k << "/" << points.at(k).getIndex() << endl;
                assert(points.at(k).getIndex()==j);
            }
        }
    }

    return P;
}

/*
 * Computes Euclidean distance ||x-y||
 */
double RBFSpline::dist(std::vector<double> x, std::vector<double> y) const
{
    assert(x.size() == y.size());
    double sum = 0.0;
    for(unsigned int i=0; i<x.size(); i++)
        sum += (x.at(i)-y.at(i))*(x.at(i)-y.at(i));
    return std::sqrt(sum);
}

/*
 * Computes Euclidean distance ||x-y||
 */
double RBFSpline::dist(DataSample x, DataSample y) const
{
    return dist(x.getX(), y.getX());
}

bool RBFSpline::dist_sort(DataSample x, DataSample y) const
{
    std::vector<double> zeros(x.getDimX(), 0);
    DataSample origin(zeros, 0.0);
    double x_dist = dist(x, origin);
    double y_dist = dist(y, origin);
    return (x_dist<y_dist);
}

} // namespace MultivariateSplines
