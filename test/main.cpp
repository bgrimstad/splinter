/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <vector>
#include <datatable.h>
#include <sstream>
#include <cmath> // abs, nextafter
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>

#include "bspline.h"
#include "pspline.h"
#include "rbfspline.h"

using std::cout;
using std::endl;

using namespace MultivariateSplines;

// Checks if a is within margin of b
bool equalsWithinRange(double a, double b, double margin = 0.0)
{
    return b - margin <= a && a <= b + margin;
}

bool is_identical(DataTable &a, DataTable &b)
{
    if(a.getNumVariables() != b.getNumVariables())
        return false;

    auto ait = a.cbegin(), bit = b.cbegin();
    for(; ait != a.cend() && bit != b.cend(); ait++, bit++)
    {
        for(unsigned int i = 0; i < a.getNumVariables(); i++)
        {
//            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getX().at(i) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getX().at(i) << " ";
            if(!equalsWithinRange(ait->getX().at(i), bit->getX().at(i)))
                return false;
        }

//            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getY().at(j) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getY().at(j) << " ";
        if(!equalsWithinRange(ait->getY(), bit->getY()))
            return false;
//        std::cout << std::endl;
    }

//    std::cout << "Finished comparing samples..." << std::endl;

    return ait == a.cend() && bit == b.cend();
}

bool test1()
{
    DataTable table;

    auto x = std::vector<double>(1);
    double y;
    for(double i = -0.3; i <= 0.3; i += 0.04)
    {
        x.at(0) = i;
        y = 2 * i;
        table.addSample(x, y);
    }

    table.save("test1.datatable");

    DataTable loadedTable;
    loadedTable.load("test1.datatable");

    remove("test1.datatable");

    return is_identical(table, loadedTable);
}

bool test2()
{
    DataTable table;

    auto x = std::vector<double>(2);
    double y;
    for(double i = -0.3; i <= 0.3; i += 0.04)
    {
        for(double j = -0.4; j <= 1.0; j += 0.08)
        {
            x.at(0) = i;
            x.at(1) = j;
            y = i * j;
            table.addSample(x, y);
        }
    }

    table.save("test2.datatable");

    DataTable loadedTable;
    loadedTable.load("test2.datatable");

    remove("test2.datatable");

    return is_identical(table, loadedTable);
}

bool test3()
{
    DataTable table;

    auto x = std::vector<double>(2);
    double y;
    for(double i = -0.3; i <= 0.3; i += 0.04)
    {
        for(double j = -0.4; j <= 1.0; j += 0.03)
        {
            x.at(0) = i;
            x.at(1) = j;
            y = i * j;
            table.addSample(x, y);
        }
    }

    table.save("test3.datatable");

    DataTable loadedTable;
    loadedTable.load("test3.datatable");

    remove("test3.datatable");

    return is_identical(table, loadedTable);
}

bool test4()
{
    DataTable table;

    auto x = std::vector<double>(4);
    double y;
    int j = 0;
    for(double i = std::numeric_limits<double>::lowest(), k = std::numeric_limits<double>::max();
        j < 10000;
        i = nextafter(i, std::numeric_limits<double>::max()), k = nextafter(k, std::numeric_limits<double>::lowest()))
    {
        x.at(0) = i;
        y = k;
        table.addSample(x, y);
        j++;
    }

    table.save("test4.datatable");

    DataTable loadedTable;
    loadedTable.load("test4.datatable");

    remove("test4.datatable");

    return is_identical(table, loadedTable);
}

bool test5()
{
    DataTable table;

    auto x = std::vector<double>(3);
    double y;
    for(double i = -0.0001; i <= 0.0001; i += 0.000001)
    {
        for(double j = -0.01; j <= 0.01; j += 0.001)
        {
            for(double k = -0.01; k <= 0.01; k += 0.001)
            {
                x.at(0) = i;
                x.at(1) = j;
                x.at(2) = k;
                y = i * j;
                table.addSample(x, y);
            }
        }
    }

    table.save("test5.datatable");

    DataTable loadedTable;
    loadedTable.load("test5.datatable");

    remove("test5.datatable");

    return is_identical(table, loadedTable);
}

bool test6()
{
    DataTable table;

    auto x = std::vector<double>(4);
    double y;
    for(double i = -0.0001; i <= 0.0001; i += 0.000001)
    {
        for(double j = -0.01; j <= 0.01; j += 0.001)
        {
            for(double k = -0.01; k <= 0.01; k += 0.001)
            {
                for(double l = -100000.0; l < 0.0; l += 13720.0)
                {
                    x.at(0) = i;
                    x.at(1) = j;
                    x.at(2) = k;
                    x.at(3) = l;
                    y = i * j;
                    table.addSample(x, y);
                }
            }
        }
    }

    table.save("test6.datatable");

    DataTable loadedTable;
    loadedTable.load("test6.datatable");

    remove("test6.datatable");

    return is_identical(table, loadedTable);
}

std::vector<double> linspace(double start, double stop, unsigned int points)
{
    std::vector<double> ret;
    double dx = 0;
    if(points > 1)
        dx = (stop - start)/(points-1);
    for(unsigned int i = 0; i < points; ++i)
        ret.push_back(start + i*dx);
    return ret;
}

// Six-hump camelback function
double f(DenseVector x)
{
    assert(x.rows() == 2);
    return (4 - 2.1*x(0)*x(0) + (1/3.)*x(0)*x(0)*x(0)*x(0))*x(0)*x(0) + x(0)*x(1) + (-4 + 4*x(1)*x(1))*x(1)*x(1);
}

void runExample()
{
    // Create new DataTable to manage samples
    DataTable samples;

    // Sample function
    auto x0_vec = linspace(0, 2, 20);
    auto x1_vec = linspace(0, 2, 20);
    DenseVector x(2);
    double y;

    for(auto x0 : x0_vec)
    {
        for(auto x1 : x1_vec)
        {
            // Sample function at x
            x(0) = x0;
            x(1) = x1;
            y = f(x);

            // Store sample
            samples.addSample(x,y);
        }
    }

    // Build B-splines that interpolate the samples
    BSpline bspline1(samples, BSplineType::LINEAR);
    BSpline bspline2(samples, BSplineType::QUADRATIC_FREE);
    BSpline bspline3(samples, BSplineType::CUBIC_FREE);

    // Build penalized B-spline (P-spline) that smooths the samples
    PSpline pspline(samples, 0.03);

    // Build radial basis function spline that interpolate the samples
    RBFSpline rbfspline(samples, RadialBasisFunctionType::THIN_PLATE_SPLINE);

    // Evaluate the splines at x = (0,0)
    x(0) = 0; x(1) = 0;

    cout << endl << endl;
    cout << "Evaluating splines at grid point x = [0,0]"        << endl;
    cout << "-------------------------------------------"       << endl;
    cout << "Function y(x):         "   << f(x)                 << endl;
    cout << "Linear B-spline:       "   << bspline1.eval(x)     << endl;
    cout << "Quadratic B-spline:    "   << bspline2.eval(x)     << endl;
    cout << "Cubic B-spline:        "   << bspline3.eval(x)     << endl;
    cout << "P-spline:              "   << pspline.eval(x)      << endl;
    cout << "Thin-plate spline:     "   << rbfspline.eval(x)    << endl;
    cout << "-------------------------------------------"       << endl;

    // Evaluate the splines at x = (1,1)
    x(0) = 1; x(1) = 1;

    cout << endl << endl;
    cout << "Evaluating splines at x = [1,1]"                   << endl;
    cout << "-------------------------------------------"       << endl;
    cout << "Function y(x):         "   << f(x)                 << endl;
    cout << "Linear B-spline:       "   << bspline1.eval(x)     << endl;
    cout << "Quadratic B-spline:    "   << bspline2.eval(x)     << endl;
    cout << "Cubic B-spline:        "   << bspline3.eval(x)     << endl;
    cout << "P-spline:              "   << pspline.eval(x)      << endl;
    cout << "Thin-plate spline:     "   << rbfspline.eval(x)    << endl;
    cout << "-------------------------------------------"       << endl;

    // Evaluate error norm
    auto x0_vec_2 = linspace(0, 2, 200);
    auto x1_vec_2 = linspace(0, 2, 200);
    std::vector<double> e_max(5, 0.0);

    for(auto x0 : x0_vec_2)
    {
        for(auto x1 : x1_vec_2)
        {
            // Sample function at x
            x(0) = x0;
            x(1) = x1;
            y = f(x);

            e_max.at(0) = std::max(e_max.at(0), std::abs(bspline1.eval(x) - y));
            e_max.at(1) = std::max(e_max.at(1), std::abs(bspline2.eval(x) - y));
            e_max.at(2) = std::max(e_max.at(2), std::abs(bspline3.eval(x) - y));
            e_max.at(3) = std::max(e_max.at(3), std::abs(pspline.eval(x) - y));
            e_max.at(4) = std::max(e_max.at(4), std::abs(rbfspline.eval(x) - y));
        }
    }

    cout << endl << endl;
    cout << "Evaluating spline errors (using max norm)  "   << endl;
    cout << "-------------------------------------------"   << endl;
    cout << "Linear B-spline:      "   << e_max.at(0)       << endl;
    cout << "Quadratic B-spline:   "   << e_max.at(1)       << endl;
    cout << "Cubic B-spline:       "   << e_max.at(2)       << endl;
    cout << "P-spline:             "   << e_max.at(3)       << endl;
    cout << "Thin-plate spline:    "   << e_max.at(4)       << endl;
    cout << "-------------------------------------------"   << endl;
}

bool compareBSplines(BSpline &bs, const BSpline &bs_orig)
{
    auto lb = bs.getDomainLowerBound();
    auto ub = bs.getDomainUpperBound();

    auto x0_vec = linspace(lb.at(0), ub.at(0), 10);
    auto x1_vec = linspace(lb.at(1), ub.at(1), 10);

    DenseVector x(2);
    for(auto x0 : x0_vec)
    {
        for(auto x1 : x1_vec)
        {
            x(0) = x0;
            x(1) = x1;

            double yb = bs.eval(x);
            double yb_orig = bs_orig.eval(x);
            if(std::abs(yb-yb_orig) > 1e-8)
            {
                cout << yb << endl;
                cout << yb_orig << endl;
                return false;
            }
        }
    }

    return true;
}

bool domainReductionTest(BSpline &bs, const BSpline &bs_orig)
{
    if(bs.getNumVariables() != 2 || bs_orig.getNumVariables() != 2)
        return false;

    // Check for error
    if(!compareBSplines(bs, bs_orig))
        return false;

    auto lb = bs.getDomainLowerBound();
    auto ub = bs.getDomainUpperBound();

    bool flag = false;
    unsigned int index = 0;
    for(; index < lb.size(); index++)
    {
        if(ub.at(index)-lb.at(index) > 1e-1)
        {
            flag = true;
            break;
        }
    }

    if(flag)
    {
        auto split = (ub.at(index) + lb.at(index))/2;

        auto lb2 = lb;
        auto ub2 = ub; ub2.at(index) = split;
        BSpline bs2(bs);
        bs2.reduceDomain(lb2, ub2);

        auto lb3 = lb; lb3.at(index) = split;
        auto ub3 = ub;
        BSpline bs3(bs);
        bs3.reduceDomain(lb3, ub3);

        return (domainReductionTest(bs2,bs_orig) && domainReductionTest(bs3,bs_orig));
    }

    return true;
}

void runRecursiveDomainReductionTest()
{
    cout << endl << endl;
    cout << "Starting recursive domain reduction test..." << endl;

    // Create new DataTable to manage samples
    DataTable samples;

    // Sample function
    auto x0_vec = linspace(0,2,20);
    auto x1_vec = linspace(0,2,20);
    DenseVector x(2);
    double y;

    for(auto x0 : x0_vec)
    {
        for(auto x1 : x1_vec)
        {
            // Sample function at x
            x(0) = x0;
            x(1) = x1;
            y = f(x);

            // Store sample
            samples.addSample(x,y);
        }
    }

    // Build B-splines that interpolate the samples
//    BSpline bspline(samples, BSplineType::LINEAR);
//    BSpline bspline(samples, BSplineType::QUADRATIC_FREE);
    BSpline bspline(samples, BSplineType::CUBIC_FREE);

    if(domainReductionTest(bspline,bspline))
        cout << "Test finished successfully!" << endl;
    else
        cout << "Test failed!" << endl;
}

void testSplineDerivative()
{
    cout << endl << endl;
    cout << "Testing spline derivative..." << endl;

    // Create new DataTable to manage samples
    DataTable samples;

    // Sample function
    double x0_lb = 0;
    double x0_ub = 2;
    double x1_lb = 0;
    double x1_ub = 2;

    auto x0_vec = linspace(x0_lb, x0_ub, 20);
    auto x1_vec = linspace(x1_lb, x1_ub, 20);
    DenseVector x(2);
    double y;

    for(auto x0 : x0_vec)
    {
        for(auto x1 : x1_vec)
        {
            // Sample function at x
            x(0) = x0;
            x(1) = x1;
            y = f(x);

            // Store sample
            samples.addSample(x,y);
        }
    }

    // Build spline that interpolate the samples
//    BSpline spline(samples, BSplineType::LINEAR);
//    BSpline spline(samples, BSplineType::QUADRATIC_FREE);
    BSpline spline(samples, BSplineType::CUBIC_FREE);
//    RBFSpline spline(samples, RadialBasisFunctionType::THIN_PLATE_SPLINE);
//    RBFSpline spline(samples, RadialBasisFunctionType::MULTIQUADRIC);

    auto x0_vec_2 = linspace(x0_lb, x0_ub, 200);
    auto x1_vec_2 = linspace(x1_lb, x1_ub, 200);

    double tol = 1e-4;  // Absolute error tolerance
    double h = 1e-8;    // Finite difference step length
    double x0_diff, x1_diff;

    for(auto x0 : x0_vec_2)
    {
        for(auto x1 : x1_vec_2)
        {
            x(0) = x0;
            x(1) = x1;

            DenseMatrix dfdx = spline.evalJacobian(x);
            if(dfdx.cols() != 2)
            {
                cout << "Test failed - check Jacobian size!" << endl;
                return;
            }

            // Finite difference in x0
            if(x(0) == x0_ub)
            {
                // Backward diff
                DenseVector dx1f = x;
                DenseVector dx1b = x; dx1b(0) = x(0)-h;
                x0_diff = (spline.eval(dx1f) - spline.eval(dx1b))/h;
            }
            else if(x(0) == x0_lb)
            {
                // Forward diff
                DenseVector dx1f = x; dx1f(0) = x(0)+h;
                DenseVector dx1b = x;
                x0_diff = (spline.eval(dx1f) - spline.eval(dx1b))/h;
            }
            else
            {
                // Central diff
                DenseVector dx1f = x; dx1f(0) = x(0)+h/2;
                DenseVector dx1b = x; dx1b(0) = x(0)-h/2;
                x0_diff = (spline.eval(dx1f) - spline.eval(dx1b))/h;
            }

            // Finite difference in x1
            if(x(1) == x1_ub)
            {
                // Backward diff
                DenseVector dx2f = x;
                DenseVector dx2b = x; dx2b(1) = x(1)-h;
                x1_diff = (spline.eval(dx2f) - spline.eval(dx2b))/h;
            }
            else if(x(1) == x1_lb)
            {
                // Forward diff
                DenseVector dx2f = x; dx2f(1) = x(1)+h;
                DenseVector dx2b = x;
                x1_diff = (spline.eval(dx2f) - spline.eval(dx2b))/h;
            }
            else
            {
                // Central diff
                DenseVector dx2f = x; dx2f(1) = x(1)+h/2;
                DenseVector dx2b = x; dx2b(1) = x(1)-h/2;
                x1_diff = (spline.eval(dx2f) - spline.eval(dx2b))/h;
            }

            if(std::abs(dfdx(0) - x0_diff) > tol
                || std::abs(dfdx(1) - x1_diff) > tol)
            {
                cout << x0 << ", " << x1 << endl;
                cout << dfdx(0) - x0_diff << endl;
                cout << dfdx(1) - x1_diff << endl;
                cout << "Test failed - check Jacobian!" << endl;
                return;
            }
        }
    }

    cout << "Test finished successfully!" << endl;
}

void run_tests()
{
    runExample();

    testSplineDerivative();

    runRecursiveDomainReductionTest();

    cout << endl << endl;
    cout << "Testing load and save functionality:       "   << endl;
    cout << "-------------------------------------------"   << endl;
    cout << "test1(): " << (test1() ? "success" : "fail")   << endl;
    cout << "test2(): " << (test2() ? "success" : "fail")   << endl;
    cout << "test3(): " << (test3() ? "success" : "fail")   << endl;
    cout << "test4(): " << (test4() ? "success" : "fail")   << endl;
    cout << "test5(): " << (test5() ? "success" : "fail")   << endl;
    cout << "test6(): " << (test6() ? "success" : "fail")   << endl;
}

int main(int argc, char **argv)
{
    try
    {
        run_tests();
    }
    catch(MultivariateSplines::Exception& e)
    {
        cout << "MS Exception - " << e.what() << endl;
    }
    catch(std::exception& e)
    {
        cout << "std::exception - " << e.what() << endl;
    }

    return 0;
}
