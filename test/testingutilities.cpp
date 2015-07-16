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
#include <test_functions.h>

using namespace std;

namespace SPLINTER
{

std::vector<std::vector<TestFunction *>> testFunctions = std::vector<std::vector<TestFunction *>>();

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
    const double one_norm_epsilon = 0.1;
    const double two_norm_epsilon = 0.1;
    const double inf_norm_epsilon = 0.1;

    return compareFunctions(exact, approx, points, one_norm_epsilon, two_norm_epsilon, inf_norm_epsilon);
}

bool compareFunctions(const Function &exact, const Function &approx, const std::vector<std::vector<double>> &points, double one_norm_epsilon, double two_norm_epsilon, double inf_norm_epsilon)
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

        /*SECTION("Function approximates the value within tolerance")*/
        {
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

            normOneVal(i) = getOneNorm(error);
            normTwoVal(i) = getTwoNorm(error);
            normInfVal(i) = getInfNorm(error);

//            REQUIRE(oneNorm(error) <= one_norm_epsilon);
//            REQUIRE(twoNorm(error) <= two_norm_epsilon);
//            REQUIRE(maxNorm(error) <= inf_norm_epsilon);
        }


        /*SECTION("Function approximates the Jacobian within tolerance")*/
        {
            auto exactJacobian = exact.evalJacobian(x);
            auto approxJacobian = approx.evalJacobian(x);
            auto errorJacobian = exactJacobian - approxJacobian;

            normOneJac(i) = getOneNorm(errorJacobian);
            normTwoJac(i) = getTwoNorm(errorJacobian);
            normInfJac(i) = getInfNorm(errorJacobian);
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


        /*SECTION("Function approximates the Hessian within tolerance")*/
        {
            auto exactHessian = exact.evalHessian(x);
            auto approxHessian = approx.evalHessian(x);
            auto errorHessian = exactHessian - approxHessian;

            normOneHes(i) = getOneNorm(errorHessian);
            normTwoHes(i) = getTwoNorm(errorHessian);
            normInfHes(i) = getInfNorm(errorHessian);

//            INFO("x: ");
//            INFO(x);
//            INFO("Exact Hessian:");
//            INFO(exactHessian);
//            INFO("Approximated Hessian:");
//            INFO(approxHessian);
//            INFO("Exact - Approx: ");
//            INFO(errorHessian);

//            CHECK(getOneNorm(errorHessian) <= one_norm_epsilon);
//            CHECK(getTwoNorm(errorHessian) <= two_norm_epsilon);
//            CHECK(getInfNorm(errorHessian) <= inf_norm_epsilon);
        }

        i++;
    }

    CHECK(getOneNorm(normOneVal) <= one_norm_epsilon);
    CHECK(getTwoNorm(normTwoVal) <= two_norm_epsilon);
    CHECK(getInfNorm(normInfVal) <= inf_norm_epsilon);

    CHECK(getOneNorm(normOneJac) <= one_norm_epsilon);
    CHECK(getTwoNorm(normTwoJac) <= two_norm_epsilon);
    CHECK(getInfNorm(normInfJac) <= inf_norm_epsilon);

    CHECK(getOneNorm(normOneHes) <= one_norm_epsilon);
    CHECK(getTwoNorm(normTwoHes) <= two_norm_epsilon);
    CHECK(getInfNorm(normInfHes) <= inf_norm_epsilon);

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


DataTable sample(const Function &func, std::vector<std::vector<double>> &points) {
    return sample(&func, points);
}

DataTable sample(const Function *func, std::vector<std::vector<double>> &points)
{
    DataTable table;

    for(auto &point : points) {
        DenseVector x = vecToDense(point);
        table.addSample(point, func->eval(x));
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

double getOneNorm(const DenseMatrix &m)
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

double getTwoNorm(const DenseMatrix &m)
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

double getInfNorm(const DenseMatrix &m)
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

#ifndef NDEBUG
    if(numPoints > 10000) {
        std::cout << "Warning: Enumerating " << numPoints << " points." << std::endl;
    }
#endif // ifndef NDEBUG

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
    auto start = std::vector<double>(dim);
    auto end = std::vector<double>(dim);
    auto numPoints = std::vector<unsigned int>(dim);


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



void setupTestFunctions() {
    // f_x_y: function of x variables and y degrees

    Var x(0, "x");
    Var y(1, "y");
    Var z(2, "z");
    Var x1(3, "x1");
    Var y1(4, "y1");
    Var z1(5, "z1");
    Var x2(6, "x2");
    Var y2(7, "y2");
    Var z2(8, "z2");

    // Functions of one variable
    auto f_1_0 = -13.37 + 0*x;
    auto f_1_1 = -5.1*x + 13.37;
    auto f_1_2 = 8.1*(x^2) - 0.2*x + 2313.1;
    auto f_1_3 = -4.5*(x^3) + 2.2*(x^2);
    auto f_1_4 = x*(4.5*(x^3) + 3*(x^2) - x);

    // Functions of two variables
    auto f_2_0 = 0.1 + 0*x*y;
    auto f_2_1 = - 5.1*x + 13.37*y;
    auto f_2_2 = 8.1*(x^2)*(y^2) - 0.2*x*y + 13.37;
    auto f_2_3 = -4.5*(x^3) + 2.2*(x^2) - (y^2) + 3;
    auto f_2_4 = x*(4.5*(x^4) - (x^2) + 3*x*y);
    auto f_2_5 = -57*(y^5)*(x^3) - 0.1*(x^4)*(y^2) + 1.1*y*(x^2) + y - 1e10;
    // Six-hump camelback function
    auto f_2_6 = (4 - 2.1*(x^2) + (1/3.)*(x^4))*(x^2) + x*y + (-4 + 4*(y^2))*(y^2) + 1.553e-1;

    // Functions of three variables
    auto f_3_0 = 6534460297 + 0*x*y*z;
    auto f_3_1 = x+y-z-1;
    auto f_3_2 = (x^2)*y*z + 3.9*(z^2) + z*(y^2) + 13.1*(z^2) - x - 10;
    auto f_3_3 = f_3_2 * f_3_1;

    // Non-polynomial (aka. nasty) functions
    auto f_nasty_x = E(z) * ((1.3^x) + (y^x));


    // First vector is the vector of nasty functions
    testFunctions.push_back(std::vector<TestFunction *>());
    testFunctions.at(0).push_back(new TestFunction(f_nasty_x));

    testFunctions.push_back(std::vector<TestFunction *>());
    testFunctions.at(1).push_back(new TestFunction(f_1_0));
    testFunctions.at(1).push_back(new TestFunction(f_1_1));
    testFunctions.at(1).push_back(new TestFunction(f_1_2));
    testFunctions.at(1).push_back(new TestFunction(f_1_3));
    testFunctions.at(1).push_back(new TestFunction(f_1_4));

    testFunctions.push_back(std::vector<TestFunction *>());
    testFunctions.at(2).push_back(new TestFunction(f_2_0));
    testFunctions.at(2).push_back(new TestFunction(f_2_1));
    testFunctions.at(2).push_back(new TestFunction(f_2_2));
    testFunctions.at(2).push_back(new TestFunction(f_2_3));
    testFunctions.at(2).push_back(new TestFunction(f_2_4));
    testFunctions.at(2).push_back(new TestFunction(f_2_5));
    testFunctions.at(2).push_back(new TestFunction(f_2_6));

    testFunctions.push_back(std::vector<TestFunction *>());
    testFunctions.at(3).push_back(new TestFunction(f_3_0));
    testFunctions.at(3).push_back(new TestFunction(f_3_1));
    testFunctions.at(3).push_back(new TestFunction(f_3_2));
    testFunctions.at(3).push_back(new TestFunction(f_3_3));
}

TestFunction *getTestFunction(int numVariables, int degree)
{
    return testFunctions.at(numVariables).at(degree);
}

std::vector<TestFunction *> getTestFunctionsOfDegree(int degree)
{
    auto testFuncs = std::vector<TestFunction *>();
    for(int i = 1; i < testFunctions.size(); ++i) {
        if(degree < testFunctions.at(i).size()) {
            testFuncs.push_back(testFunctions.at(i).at(degree));
        }
    }
    return testFuncs;
}

std::vector<TestFunction *> getTestFunctionWithNumVariables(int numVariables)
{
    return testFunctions.at(numVariables);
}

std::vector<TestFunction *> getNiceTestFunctions()
{
    auto testFuncs = std::vector<TestFunction *>();
    for(int i = 1; i < testFunctions.size(); ++i) {
        for(int j = 0; j < testFunctions.at(i).size(); ++j) {
            testFuncs.push_back(testFunctions.at(i).at(j));
        }
    }
    return testFuncs;
}

std::vector<TestFunction *> getNastyTestFunctions()
{
    return testFunctions.at(0);
}


} // namespace SPLINTER
