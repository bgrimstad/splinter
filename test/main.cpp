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
    double x_l = 0;
    double x_u = 2;
    int num_samples = 20;
    double delta_x = (x_u - x_l)/(num_samples-1);
    DenseVector x(2);
    double y;

    for(int i = 0; i < num_samples; i++)
    {
        for(int j = 0; j < num_samples; j++)
        {
            // Sample function at x
            x(0) = x_l + i*delta_x;
            x(1) = x_l + j*delta_x;
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
    std::vector<double> e_max(5,0.0);
    int num_samples_2 = 200;
    double delta_x_2 = (x_u - x_l)/(num_samples_2-1);

    for(int i = 0; i < num_samples_2; i++)
    {
        for(int j = 0; j < num_samples_2; j++)
        {
            // Sample function at x
            x(0) = i*delta_x_2;
            x(1) = j*delta_x_2;
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
    cout << endl << endl;

}

int main(int argc, char **argv)
{
    runExample();

    cout << "test1(): " << (test1() ? "success" : "fail") << endl;
    cout << "test2(): " << (test2() ? "success" : "fail") << endl;
    cout << "test3(): " << (test3() ? "success" : "fail") << endl;
    cout << "test4(): " << (test4() ? "success" : "fail") << endl;
    cout << "test5(): " << (test5() ? "success" : "fail") << endl;
    cout << "test6(): " << (test6() ? "success" : "fail") << endl;

    return 0;
}
