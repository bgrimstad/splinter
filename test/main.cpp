/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#define CATCH_CONFIG_RUNNER
#include <Catch.h>
#include <testfunction.h>
#include <testfunctions.h>

int main(int argc, char *argv[])
{
    setupTestFunctions();

    int result = Catch::Session().run(argc, argv);

    tearDownTestFunctions();

    return result;
}

//
//double kroneckerTestFunction(DenseVector x)
//{
////    assert(x.rows() == dim);
//
//    double y = 1 + (.1 + 0.5*x(0) - x(0)*x(0) - 0.33*x(0)*x(0)*x(0))*(.1 + 0.5*x(1) - x(1)*x(1) - 2*x(1)*x(1)*x(1))*(.1 - 0.5*x(2) + x(2)*x(2) + 2*x(2)*x(2)*x(2))*(.1 - 0.5*x(3) + x(3)*x(3) + 2*x(3)*x(3)*x(3));
//    return y;
//}
//
///*
// * This example was used to test the speed difference
// * of Eigens Kronecker product and the custom-made
// * myKroneckerProduct() when refining a B-spline.
// * Note that the Kronecker product function call
// * have to be manually switched between test runs.
// */
//void kroneckerTest()
//{
//    std::vector<double> lb = {-1,-1,-1,-1};
//    std::vector<double> ub = {2,2,2,2};
//
//    Sample table;
//    int n = 6;
//    for (int i0 = 0; i0 < n; i0++)
//    {
//        for (int i1 = 0; i1 < n; i1++)
//        {
//            for (int i2 = 0; i2 < n; i2++)
//            {
//                for (int i3 = 0; i3 < n; i3++)
//                {
//                    double x0 = lb.at(0) + i0*(ub.at(0)-lb.at(0))/(n-1);
//                    double x1 = lb.at(1) + i1*(ub.at(1)-lb.at(1))/(n-1);
//                    double x2 = lb.at(2) + i2*(ub.at(2)-lb.at(2))/(n-1);
//                    double x3 = lb.at(3) + i3*(ub.at(3)-lb.at(3))/(n-1);
//
//                    DenseVector xv(4);
//                    xv(0) = x0;
//                    xv(1) = x1;
//                    xv(2) = x2;
//                    xv(3) = x3;
//
//                    double y = kroneckerTestFunction(xv);
//                    table.addSample(xv,y);
//                }
//            }
//        }
//    }
//
//    cout << "Creating large 4-D B-spline" << endl;
//    BSpline bs(table, BSplineType::CUBIC);
//    BSpline bs2(bs);
//
//    //Timer timer;
//
//    std::vector<double> lb2 = {-.5,-.5,-.5,-.5};
//    std::vector<double> ub2 = {1,1,1,1};
//
//    cout << "Reducing domain" << endl;
//    //timer.start();
//    bs.reduceDomain(lb2,ub2);
//    //timer.stop();
//    // Old insertion method: 538, 536
//    // New insertion method: > 25000!
//    //cout << "Time (ms): " << timer.getMilliSeconds() << endl;
//
//    cout << "Doing error check!" << endl;
//    for (double i0 = lb2.at(0); i0 <= ub2.at(0); i0 += 0.1)
//    {
//        for (double i1 = lb2.at(1); i1 <= ub2.at(1); i1 += 0.1)
//        {
//            for (double i2 = lb2.at(2); i2 <= ub2.at(2); i2 += 0.1)
//            {
//                for (double i3 = lb2.at(3); i3 <= ub2.at(3); i3 += 0.1)
//                {
//                    DenseVector x(4);
//                    x(0) = i0;
//                    x(1) = i1;
//                    x(2) = i2;
//                    x(3) = i3;
//
//                    auto y1 = bs.eval(x);
//                    auto y2 = bs2.eval(x);
//
//                    auto dy = y1-y2;
//                    //dymax(0) += 1e-12;
//
//                    double eps = std::numeric_limits<double>::epsilon();
//                    eps = 1e-14;
//                    if (dy > eps || dy < -eps)
//                    {
//                        cout << "Error detected: " << dy << endl;
//                        exit(1);
//                    }
//
//                }
//            }
//        }
//    }
//    cout << "Error check ran successfully!" << endl;
//
//    // Check knot vectors after domain reduction
////    std::vector< std::vector<double> > knots2 = bs2.getKnotVectors();
////    std::vector< std::vector<double> > knots = bs.getKnotVectors();
//
////    cout << "Original knots: " << endl;
////    cout << "---------------------" << endl;
////    printVector(knots2);
//
////    cout << "New knots: " << endl;
////    cout << "---------------------" << endl;
////    printVector(knots);
//
//}
//
//void localRefinementTest()
//{
////    DenseMatrix coeffs = DenseMatrix::Ones(1,4);
////    std::vector<std::vector<double>> knots = {{1,1,2,3,4,4}};
//
//    DenseMatrix coeffs = DenseMatrix::Ones(1,2);
//    std::vector<std::vector<double>> knots = {{1,1,1.000000001,1.000000001}};
//
//    std::vector<unsigned int> degs = {1};
//    BSpline bs(coeffs, knots, degs);
//
//    /*
//     * x = 1 => knot inserted at 1.5
//     * x = 2 => knot inserted at 1.5
//     * x = 2.0001 => knot inserted at 2.5
//     * x = 2.01 => knot inserted at 2.001
//     * x = 4 => knot inserted at 3.5
//     */
//    DenseVector x(1); x(0) = 1;
//    bs.localKnotRefinement(x);
//
//    auto knots2 = bs.getKnotVectors();
//
//    for (auto k : knots2.at(0))
//        cout << k << ", ";
//}
//
