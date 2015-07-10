/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "testingutilities.h"
#include <iostream>
#include <Catch.h>

using namespace std;

namespace SPLINTER
{

// Checks if a is within margin of b
bool equalsWithinRange(double a, double b, double margin)
{
    return b - margin <= a && a <= b + margin;
}

bool compareFunctions(const Function &f1, const Function &f2)
{
    int dim = f1.getNumVariables();
    auto start = std::vector<double>(dim);
    auto end = std::vector<double>(dim);
    auto numPoints = std::vector<unsigned int>(dim);

    // Try to avoid testing at "nice" values,
    // as that is likely where the function was sampled.
    for(int i = 0; i < dim; i++) {
        start.at(i) = -10.0;
        end.at(i) = 10.0;
        numPoints.at(i) = 10;
    }

    auto points = linspace(start, end, numPoints);

    return compareFunctions(f1, f2, points);
}

bool compareFunctions(const Function &exact, const Function &approx, const std::vector<std::vector<double>> &points)
{
    // Max error in function value
    const double val_epsilon = 0.01;
    const double one_norm_epsilon = 0.1;
    const double two_norm_epsilon = 0.1;
    const double inf_norm_epsilon = 0.1;

    return compareFunctions(exact, approx, points, val_epsilon, one_norm_epsilon, two_norm_epsilon, inf_norm_epsilon);
}
bool compareFunctions(const Function &exact, const Function &approx, const std::vector<std::vector<double>> &points, double val_epsilon, double one_norm_epsilon, double two_norm_epsilon, double inf_norm_epsilon)
{
    REQUIRE(exact.getNumVariables() == approx.getNumVariables());

    DenseVector normOneVal(points.size());
    DenseVector normTwoVal(points.size());
    DenseVector normInfVal(points.size());

    DenseVector normOneJac(points.size());
    DenseVector normTwoJac(points.size());
    DenseVector normInfJac(points.size());

    DenseVector normOneHes(points.size());
    DenseVector normTwoHes(points.size());
    DenseVector normInfHes(points.size());

    int i = 0;
    for (auto &point : points) {
        DenseVector x = vecToDense(point);

//        INFO("Evaluation point: " << pretty_print(x));

        /*SECTION("Function approximates the value within tolerance")*/ {
            DenseMatrix exactValue(1,1);
            exactValue(0,0) = exact.eval(x);
            DenseMatrix approxValue(1,1);
            approxValue(0,0) = approx.eval(x);
            DenseMatrix error = exactValue - approxValue;

//            INFO("Exact value:");
//            INFO(exactValue);
//            INFO("Approximated value:");
//            INFO(approxValue);
//            INFO("Exact - approx:");
//            INFO(error);

            normOneVal(i) = oneNorm(error);
            normTwoVal(i) = twoNorm(error);
            normInfVal(i) = maxNorm(error);

//            REQUIRE(oneNorm(error) <= one_norm_epsilon);
//            REQUIRE(twoNorm(error) <= two_norm_epsilon);
//            REQUIRE(maxNorm(error) <= inf_norm_epsilon);
        }


        /*SECTION("Function approximates the Jacobian within tolerance")*/ {
            auto exactJacobian = exact.evalJacobian(x);
            auto approxJacobian = approx.evalJacobian(x);
            auto errorJacobian = exactJacobian - approxJacobian;

            normOneJac(i) = oneNorm(errorJacobian);
            normTwoJac(i) = twoNorm(errorJacobian);
            normInfJac(i) = maxNorm(errorJacobian);
//            INFO("Exact Jacobian:");
//            INFO(exactJacobian);
//            INFO("Approximated Jacobian:");
//            INFO(approxJacobian);
//            INFO("Exact - Approx: ");
//            INFO(errorJacobian);
//
//            REQUIRE(oneNorm(errorJacobian) <= one_norm_epsilon);
//            REQUIRE(twoNorm(errorJacobian) <= two_norm_epsilon);
//            REQUIRE(maxNorm(errorJacobian) <= inf_norm_epsilon);
        }


        /*SECTION("Function approximates the Hessian within tolerance")*/ {
            auto exactHessian = exact.evalHessian(x);
            auto approxHessian = approx.evalHessian(x);
            auto errorHessian = exactHessian - approxHessian;

            normOneHes(i) = oneNorm(errorHessian);
            normTwoHes(i) = twoNorm(errorHessian);
            normInfHes(i) = maxNorm(errorHessian);

//            INFO("Exact Hessian:");
//            INFO(exactHessian);
//            INFO("Approximated Hessian:");
//            INFO(approxHessian);
//            INFO("Exact - Approx: ");
//            INFO(errorHessian);
//
//            CHECK(oneNorm(errorHessian) <= one_norm_epsilon);
//            CHECK(twoNorm(errorHessian) <= two_norm_epsilon);
//            CHECK(maxNorm(errorHessian) <= inf_norm_epsilon);
        }

        i++;
    }

    CHECK(oneNorm(normOneVal) <= one_norm_epsilon);
    CHECK(oneNorm(normOneJac) <= one_norm_epsilon);
    CHECK(oneNorm(normOneHes) <= one_norm_epsilon);

    CHECK(twoNorm(normTwoVal) <= two_norm_epsilon);
    CHECK(twoNorm(normTwoJac) <= two_norm_epsilon);
    CHECK(twoNorm(normTwoHes) <= two_norm_epsilon);

    CHECK(maxNorm(normInfVal) <= inf_norm_epsilon);
    CHECK(maxNorm(normInfJac) <= inf_norm_epsilon);
    CHECK(maxNorm(normInfHes) <= inf_norm_epsilon);

    return true;
}

bool compareBSplines(const BSpline &left, const BSpline &right)
{
    auto left_lb = left.getDomainLowerBound();
    auto left_ub = left.getDomainUpperBound();
    auto right_lb = right.getDomainLowerBound();
    auto right_ub = right.getDomainUpperBound();

    REQUIRE(left_lb.size() == left_ub.size());
    REQUIRE(left_ub.size() == right_lb.size());
    REQUIRE(right_lb.size() == right_ub.size());

    int dim = left_lb.size();

    auto points = linspace(dim);

    for(int i = 0; i < dim; i++) {
        REQUIRE(left_lb.at(i) == right_lb.at(i));
        REQUIRE(left_ub.at(i) == right_ub.at(i));
    }

    return compareFunctions(left, right, points);

//    auto x0_vec = linspace(lb.at(0), ub.at(0), 10);
//    auto x1_vec = linspace(lb.at(1), ub.at(1), 10);
//
//    DenseVector x(2);
//    for (auto x0 : x0_vec)
//    {
//        for (auto x1 : x1_vec)
//        {
//            x(0) = x0;
//            x(1) = x1;
//
//            double yb = bs.eval(x);
//            double yb_orig = bs_orig.eval(x);
//            if (std::abs(yb-yb_orig) > 1e-8)
//            {
//                cout << yb << endl;
//                cout << yb_orig << endl;
//                return false;
//            }
//        }
//    }
//
//    return true;
}

bool compareDataTables(DataTable &a, DataTable &b)
{
    if (a.getNumVariables() != b.getNumVariables())
        return false;

    auto ait = a.cbegin(), bit = b.cbegin();
    for (; ait != a.cend() && bit != b.cend(); ait++, bit++)
    {
        for (unsigned int i = 0; i < a.getNumVariables(); i++)
        {
//            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getX().at(i) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getX().at(i) << " ";
            if (!equalsWithinRange(ait->getX().at(i), bit->getX().at(i)))
                return false;
        }

//            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getY().at(j) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getY().at(j) << " ";
        if (!equalsWithinRange(ait->getY(), bit->getY()))
            return false;
//        std::cout << std::endl;
    }

//    std::cout << "Finished comparing samples..." << std::endl;

    return ait == a.cend() && bit == b.cend();
}

DataTable sample(const Function &func, std::vector<std::vector<double>> &points)
{
    DataTable table;

    for(auto &point : points) {
        DenseVector x = vecToDense(point);
        table.addSample(point, func.eval(x));
    }

    return table;
}

std::vector<double> linspace(double start, double stop, unsigned int num)
{
    std::vector<double> ret;
    double dx = 0;
    if (num > 1)
        dx = (stop - start)/(num-1);
    for (unsigned int i = 0; i < num; ++i)
        ret.push_back(start + i*dx);
    return ret;
}

double sixHumpCamelBack(DenseVector x)
{
    assert(x.rows() == 2);
    return (4 - 2.1*x(0)*x(0) + (1/3.)*x(0)*x(0)*x(0)*x(0))*x(0)*x(0) + x(0)*x(1) + (-4 + 4*x(1)*x(1))*x(1)*x(1);
}

double oneNorm(const DenseMatrix &m)
{
    return m.lpNorm<1>();
    double norm = 0.0;
    for(int i = 0; i < m.rows(); i++) {
        for(int j = 0; j < m.cols(); j++) {
            norm += std::abs(m(i, j));
        }
    }

    return norm;
}

double twoNorm(const DenseMatrix &m)
{
    return m.lpNorm<2>();
    double norm = 0.0;
    for(int i = 0; i < m.rows(); i++) {
        for(int j = 0; j < m.cols(); j++) {
            norm += std::pow(m(i, j), 2);
        }
    }

    norm = std::sqrt(norm);

    return norm;
}

double maxNorm(const DenseMatrix &m)
{
    return m.lpNorm<Eigen::Infinity>();
    double norm = 0.0;
    for(int i = 0; i < m.rows(); i++) {
        for(int j = 0; j < m.cols(); j++) {
            norm = std::max(std::abs(m(i, j)), norm);
        }
    }

    return norm;
}

double log(double base, double x)
{
    return std::log(x) / std::log(base);
}

// Enumerates all permutations
std::vector<std::vector<double>> linspace(std::vector<double> start, std::vector<double> end, std::vector<unsigned int> points)
{
    assert(start.size() == end.size() && end.size() == points.size());

    unsigned int nDims = start.size();
    unsigned int numPoints = 1;
    for(int i = 0; i < nDims; i++) {
        numPoints *= points.at(i);
    }

    if(true || numPoints > 10000) {
        std::cout << "Warning: Enumerating " << numPoints << " points." << std::endl;
    }

    auto result = std::vector<std::vector<double>>(numPoints);

    auto dx = std::vector<double>(nDims);
    auto it = std::vector<unsigned int>(nDims);
    auto cur = std::vector<double>(nDims);

    for(unsigned int i = 0; i < nDims; i++) {
        cur.at(i) = start.at(i);
        it.at(i) = 0;
        dx.at(i) = (end.at(i) - start.at(i)) / (points.at(i) - 1);
    }

    unsigned int curDim = 0;
    unsigned int i = 0;
    // Add the start vector
    result.at(i++) = std::vector<double>(start);
    while(true) {
        curDim = 0;
        while(curDim < nDims && it.at(curDim)+1 >= points.at(curDim)) {
            it.at(curDim) = 0;
            cur.at(curDim) = start.at(curDim) + dx.at(curDim) * it.at(curDim);
            curDim++;
        }
        // If we're unable to find a number that can be increased,
        // it means were at the highest number
        // (for 4 digits in decimal that is 9999)
        if(curDim >= nDims) {
            break;
        }

        it.at(curDim)++;
        cur.at(curDim) = start.at(curDim) + dx.at(curDim) * it.at(curDim);

        result.at(i++) = std::vector<double>(cur);
    }

    assert(i == numPoints);

    return result;
}

std::vector<std::vector<double>> linspace(int dim, double start, double end, unsigned int points)
{
    auto startVec = std::vector<double>(dim);
    auto endVec = std::vector<double>(dim);
    auto pointsVec = std::vector<unsigned int>(dim);

    for(int i = 0; i < dim; i++) {
        startVec.at(i) = start;
        endVec.at(i) = end;
        pointsVec.at(i) = points;
    }

    return linspace(startVec, endVec, pointsVec);
}

std::vector<std::vector<double>> linspace(int dim)
{
    // Total number of points to evaluate
    // Will be distributed evenly per dimension
    const int totSamples = 10000;

    // Take the dim'th root to get number of samples per dimension
    unsigned int pointsPerDim = std::pow(totSamples, 1.0/dim);

    // Clamp so that pointsPerDim >= 2
    if(pointsPerDim < 2) {
        pointsPerDim = 2;
    }

    return linspace(dim, pointsPerDim);
}

std::vector<std::vector<double>> linspace(int dim, unsigned int pointsPerDim)
{
    auto start = std::vector<double>();
    auto end = std::vector<double>();
    auto numPoints = std::vector<unsigned int>();


    for(int i = 0; i < dim; i++) {

        start.at(i) = -10.0;
        end.at(i) = 10.0;
        numPoints.at(i) = pointsPerDim;
    }

    return linspace(start, end, numPoints);
}


std::vector<double> denseToVec(const DenseVector &dense)
{
    auto vec = std::vector<double>(dense.size());
    for(int i = 0; i < dense.size(); i++) {
        vec.at(i) = dense(i);
    }

    return vec;
}

DenseVector vecToDense(const std::vector<double> &vec)
{
    DenseVector dense(vec.size());
    for(int i = 0; i < vec.size(); i++) {
        dense(i) = vec.at(i);
    }

    return dense;
}

std::string pretty_print(const DenseVector &denseVec)
{
    std::string str("[");
    for(int i = 0; i < denseVec.rows(); i++) {
        str += to_string(denseVec(i));
        if(i + 1 < denseVec.rows()) {
            str += "; ";
        }
    }
    str += "]";

    return str;
}

} // namespace SPLINTER