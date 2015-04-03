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
#include "testingutilities.h"

using std::cout;
using std::endl;

using namespace MultivariateSplines;


void runExample()
{
    // Create new DataTable to manage samples
    DataTable samples;

    // Sample function
    auto x0_vec = linspace(0, 2, 20);
    auto x1_vec = linspace(0, 2, 20);
    DenseVector x(2);
    double y;

    for (auto x0 : x0_vec)
    {
        for (auto x1 : x1_vec)
        {
            // Sample function at x
            x(0) = x0;
            x(1) = x1;
            y = sixHumpCamelBack(x);

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
    cout << "Function y(x):         "   << sixHumpCamelBack(x)  << endl;
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
    cout << "Function y(x):         "   << sixHumpCamelBack(x)  << endl;
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

    for (auto x0 : x0_vec_2)
    {
        for (auto x1 : x1_vec_2)
        {
            // Sample function at x
            x(0) = x0;
            x(1) = x1;
            y = sixHumpCamelBack(x);

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

bool domainReductionTest(BSpline &bs, const BSpline &bs_orig)
{
    if (bs.getNumVariables() != 2 || bs_orig.getNumVariables() != 2)
        return false;

    // Check for error
    if (!compareBSplines(bs, bs_orig))
        return false;

    auto lb = bs.getDomainLowerBound();
    auto ub = bs.getDomainUpperBound();

    bool flag = false;
    unsigned int index = 0;
    for (; index < lb.size(); index++)
    {
        if (ub.at(index)-lb.at(index) > 1e-1)
        {
            flag = true;
            break;
        }
    }

    if (flag)
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

    for (auto x0 : x0_vec)
    {
        for (auto x1 : x1_vec)
        {
            // Sample function at x
            x(0) = x0;
            x(1) = x1;
            y = sixHumpCamelBack(x);

            // Store sample
            samples.addSample(x,y);
        }
    }

    // Build B-splines that interpolate the samples
//    BSpline bspline(samples, BSplineType::LINEAR);
//    BSpline bspline(samples, BSplineType::QUADRATIC_FREE);
    BSpline bspline(samples, BSplineType::CUBIC_FREE);

    if (domainReductionTest(bspline,bspline))
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

    for (auto x0 : x0_vec)
    {
        for (auto x1 : x1_vec)
        {
            // Sample function at x
            x(0) = x0;
            x(1) = x1;
            y = sixHumpCamelBack(x);

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

    for (auto x0 : x0_vec_2)
    {
        for (auto x1 : x1_vec_2)
        {
            x(0) = x0;
            x(1) = x1;

            DenseMatrix dfdx = spline.evalJacobian(x);
            if (dfdx.cols() != 2)
            {
                cout << "Test failed - check Jacobian size!" << endl;
                return;
            }

            // Finite difference in x0
            if (x(0) == x0_ub)
            {
                // Backward diff
                DenseVector dx1f = x;
                DenseVector dx1b = x; dx1b(0) = x(0)-h;
                x0_diff = (spline.eval(dx1f) - spline.eval(dx1b))/h;
            }
            else if (x(0) == x0_lb)
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
            if (x(1) == x1_ub)
            {
                // Backward diff
                DenseVector dx2f = x;
                DenseVector dx2b = x; dx2b(1) = x(1)-h;
                x1_diff = (spline.eval(dx2f) - spline.eval(dx2b))/h;
            }
            else if (x(1) == x1_lb)
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

            if (std::abs(dfdx(0) - x0_diff) > tol
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

double kroneckerTestFunction(DenseVector x)
{
    int dim = 4;
    assert(x.rows() == dim);

    double y = 1 + (.1 + 0.5*x(0) - x(0)*x(0) - 0.33*x(0)*x(0)*x(0))*(.1 + 0.5*x(1) - x(1)*x(1) - 2*x(1)*x(1)*x(1))*(.1 - 0.5*x(2) + x(2)*x(2) + 2*x(2)*x(2)*x(2))*(.1 - 0.5*x(3) + x(3)*x(3) + 2*x(3)*x(3)*x(3));
    return y;
}

/*
 * This example was used to test the speed difference
 * of Eigens Kronecker product and the custom-made
 * myKroneckerProduct() when refining a B-spline.
 * Note that the Kronecker product function call
 * have to be manually switched between test runs.
 */
void kroneckerTest()
{
    int dim = 4;

    std::vector<double> lb = {-1,-1,-1,-1};
    std::vector<double> ub = {2,2,2,2};

    DataTable table;
    int n = 6;
    for (int i0 = 0; i0 < n; i0++)
    {
        for (int i1 = 0; i1 < n; i1++)
        {
            for (int i2 = 0; i2 < n; i2++)
            {
                for (int i3 = 0; i3 < n; i3++)
                {
                    double x0 = lb.at(0) + i0*(ub.at(0)-lb.at(0))/(n-1);
                    double x1 = lb.at(1) + i1*(ub.at(1)-lb.at(1))/(n-1);
                    double x2 = lb.at(2) + i2*(ub.at(2)-lb.at(2))/(n-1);
                    double x3 = lb.at(3) + i3*(ub.at(3)-lb.at(3))/(n-1);

                    DenseVector xv(4);
                    xv(0) = x0;
                    xv(1) = x1;
                    xv(2) = x2;
                    xv(3) = x3;

                    double y = kroneckerTestFunction(xv);
                    table.addSample(xv,y);
                }
            }
        }
    }

    cout << "Creating large 4-D B-spline" << endl;
    BSpline bs(table, BSplineType::CUBIC_FREE);
    BSpline bs2(bs);

    //Timer timer;

    std::vector<double> lb2 = {-.5,-.5,-.5,-.5};
    std::vector<double> ub2 = {1,1,1,1};

    cout << "Reducing domain" << endl;
    //timer.start();
    bs.reduceDomain(lb2,ub2);
    //timer.stop();
    // Old insertion method: 538, 536
    // New insertion method: > 25000!
    //cout << "Time (ms): " << timer.getMilliSeconds() << endl;

    cout << "Doing error check!" << endl;
    for (double i0 = lb2.at(0); i0 <= ub2.at(0); i0 += 0.1)
    {
        for (double i1 = lb2.at(1); i1 <= ub2.at(1); i1 += 0.1)
        {
            for (double i2 = lb2.at(2); i2 <= ub2.at(2); i2 += 0.1)
            {
                for (double i3 = lb2.at(3); i3 <= ub2.at(3); i3 += 0.1)
                {
                    DenseVector x(4);
                    x(0) = i0;
                    x(1) = i1;
                    x(2) = i2;
                    x(3) = i3;

                    auto y1 = bs.eval(x);
                    auto y2 = bs2.eval(x);

                    auto dy = y1-y2;
                    //dymax(0) += 1e-12;

                    double eps = std::numeric_limits<double>::epsilon();
                    eps = 1e-14;
                    if (dy > eps || dy < -eps)
                    {
                        cout << "Error detected: " << dy << endl;
                        exit(1);
                    }

                }
            }
        }
    }
    cout << "Error check ran successfully!" << endl;

    // Check knot vectors after domain reduction
//    std::vector< std::vector<double> > knots2 = bs2.getKnotVectors();
//    std::vector< std::vector<double> > knots = bs.getKnotVectors();

//    cout << "Original knots: " << endl;
//    cout << "---------------------" << endl;
//    printVector(knots2);

//    cout << "New knots: " << endl;
//    cout << "---------------------" << endl;
//    printVector(knots);

}

void localRefinementTest()
{
    DenseMatrix coeffs = DenseMatrix::Ones(1,4);
    std::vector<std::vector<double>> knots = {{1,1,2,3,4,4}};
    std::vector<uint> degs = {1};
    BSpline bs(coeffs, knots, degs);

    /*
     * x = 1 => knot inserted at 1.5
     * x = 2 => knot inserted at 1.5
     * x = 2.0001 => knot inserted at 2.5
     * x = 2.01 => knot inserted at 2.001
     * x = 4 => knot inserted at 3.5
     */
    DenseVector x(1); x(0) = 2.01;
    bs.localRefinement(x);

    auto knots2 = bs.getKnotVectors();

    for (auto k : knots2.at(0))
        cout << k << ", ";
}

void run_tests()
{
    runExample();

    testSplineDerivative();

    runRecursiveDomainReductionTest();
}

int main(int argc, char **argv)
{
    try
    {
        run_tests();
    }
    catch (MultivariateSplines::Exception& e)
    {
        cout << "MS Exception - " << e.what() << endl;
    }
    catch (std::exception& e)
    {
        cout << "std::exception - " << e.what() << endl;
    }

    return 0;
}
