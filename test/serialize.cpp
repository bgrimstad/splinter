/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <datatable.h>
#include <datasample.h>
#include "testingutilities.h"
#include <iostream>
#include <bspline.h>
#include <pspline.h>
#include <radialbasisfunction.h>
#include <polynomialregression.h>
#include <serializer.h>
#include <chrono>

using namespace std;
using namespace SPLINTER;

bool serializeDataTable1()
{
    DataTable table;

    auto x = std::vector<double>(1);
    double y;
    for (double i = -0.3; i <= 0.3; i += 0.04)
    {
        x.at(0) = i;
        y = 2 * i;
        table.addSample(x, y);
    }

    table.save("test1.datatable");

    DataTable loadedTable("test1.datatable");

    remove("test1.datatable");

    return compareDataTables(table, loadedTable);
}

bool serializeDataTable2()
{
    DataTable table;

    auto x = std::vector<double>(2);
    double y;
    for (double i = -0.3; i <= 0.3; i += 0.04)
    {
        for (double j = -0.4; j <= 1.0; j += 0.08)
        {
            x.at(0) = i;
            x.at(1) = j;
            y = i * j;
            table.addSample(x, y);
        }
    }

    table.save("test2.datatable");

    DataTable loadedTable("test2.datatable");

    remove("test2.datatable");

    return compareDataTables(table, loadedTable);
}

bool serializeDataTable3()
{
    DataTable table;

    auto x = std::vector<double>(2);
    double y;
    for (double i = -0.3; i <= 0.3; i += 0.04)
    {
        for (double j = -0.4; j <= 1.0; j += 0.03)
        {
            x.at(0) = i;
            x.at(1) = j;
            y = i * j;
            table.addSample(x, y);
        }
    }

    table.save("test3.datatable");

    DataTable loadedTable("test3.datatable");

    remove("test3.datatable");

    return compareDataTables(table, loadedTable);
}

bool serializeDataTable4()
{
    DataTable table;

    auto x = std::vector<double>(4);
    double y;
    int j = 0;
    for (double i = std::numeric_limits<double>::lowest(), k = std::numeric_limits<double>::max();
        j < 10000;
        i = nextafter(i, std::numeric_limits<double>::max()), k = nextafter(k, std::numeric_limits<double>::lowest()))
    {
        x.at(0) = i;
        y = k;
        table.addSample(x, y);
        j++;
    }

    table.save("test4.datatable");

    DataTable loadedTable("test4.datatable");

    remove("test4.datatable");

    return compareDataTables(table, loadedTable);
}

bool serializeDataTable5()
{
    DataTable table;

    auto x = std::vector<double>(3);
    double y;
    for (double i = -0.0001; i <= 0.0001; i += 0.000001)
    {
        for (double j = -0.01; j <= 0.01; j += 0.001)
        {
            for (double k = -0.01; k <= 0.01; k += 0.001)
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

    DataTable loadedTable("test5.datatable");

    remove("test5.datatable");

    return compareDataTables(table, loadedTable);
}

bool serializeDataTable6()
{
    DataTable table;

    auto x = std::vector<double>(4);
    double y;
    for (double i = -0.0001; i <= 0.0001; i += 0.000001)
    {
        for (double j = -0.01; j <= 0.01; j += 0.001)
        {
            for (double k = -0.01; k <= 0.01; k += 0.001)
            {
                for (double l = -100000.0; l < -60000.0; l += 13720.0)
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

    DataTable loadedTable("test6.datatable");

    remove("test6.datatable");

    return compareDataTables(table, loadedTable);
}

bool serializeBSpline1()
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
    BSpline bspline(samples, BSplineType::LINEAR);
    bspline.save("saveTest1.bspline");
    BSpline loadedBspline("saveTest1.bspline");
    remove("saveTest1.bspline");

    return compareBSplines(bspline, loadedBspline);
}

bool serializePSpline1()
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
    // Build P-spline
    PSpline pspline(samples);
    pspline.save("saveTest1.pspline");
    PSpline loadedPspline("saveTest1.pspline");
    remove("saveTest1.pspline");

    return compareBSplines(pspline, loadedPspline)
           && pspline.getLambda() == loadedPspline.getLambda();
}

bool serializeRBF1()
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
    // Build RBF
    RadialBasisFunction rbf(samples, RadialBasisFunctionType::THIN_PLATE_SPLINE);
    rbf.save("saveTest1.rbf");
    RadialBasisFunction loadedRbf("saveTest1.rbf");
    remove("saveTest1.rbf");

    return compareFunctions(rbf, loadedRbf);
}

int main()
{
    /*cout << "Serializer test: " << std::endl;

    auto vec = std::vector<std::vector<double>>(10000);
    for(int i = 0; i < vec.size(); i++) {
        for(int j = 0; j < 1000; j++) {
            vec.at(i).push_back(j*10);
        }
    }

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    StreamType stream;
    serialize(vec, stream);
    serialize(vec, stream);
    save_to_file("test2.serialized", stream);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    std::cout << "Old method serialize = " << duration_cast<duration<double>>(t2 - t1).count() << std::endl;

    StreamType stream2 = load_from_file("test2.serialized");
    auto it = stream2.cbegin();
    auto vec3 = deserialize<std::vector<std::vector<double>>>(it, stream2.cend());
    auto vec5 = deserialize<std::vector<std::vector<double>>>(it, stream2.cend());

    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    std::cout << "Old method deserialize = " << duration_cast<duration<double>>(t3 - t2).count() << std::endl;

    Serializer s;
    s.serialize(vec);
    s.serialize(vec);
    s.saveToFile("test.serialized");

    high_resolution_clock::time_point t4 = high_resolution_clock::now();
    std::cout << "New method serialize = " << duration_cast<duration<double>>(t4 - t3).count() << std::endl;

    s.loadFromFile("test.serialized");
    std::vector<std::vector<double>> vec2;
    std::vector<std::vector<double>> vec4;
    s.deserialize(vec2);
    s.deserialize(vec4);

    high_resolution_clock::time_point t5 = high_resolution_clock::now();
    std::cout << "New method deserialize = " << duration_cast<duration<double>>(t5 - t4).count() << std::endl;


    for(int i = 0; i < vec.size(); i++) {
        for(int j = 0; j < vec.at(i).size(); j++) {
            if(vec.at(i).at(j) != vec2.at(i).at(j)) {
                std::cout << "vec not equal to vec2!" << std::endl;
            }
            if(vec.at(i).at(j) != vec3.at(i).at(j)) {
                std::cout << "vec not equal to vec3!" << std::endl;
            }
            if(vec.at(i).at(j) != vec4.at(i).at(j)) {
                std::cout << "vec not equal to vec4!" << std::endl;
            }
            if(vec.at(i).at(j) != vec5.at(i).at(j)) {
                std::cout << "vec not equal to vec5!" << std::endl;
            }
        }
    }

    exit(0);*/

    using namespace std::chrono;
    cout << endl << endl;
    cout << "Testing load and save functionality:                                   "   << endl;
    cout << "-----------------------------------------------------------------------"   << endl;
    cout << "DataTables:                                                            "   << endl;
    cout << "serializeDataTable1(): " << (serializeDataTable1() ? "success" : "fail")   << endl;
    cout << "serializeDataTable2(): " << (serializeDataTable2() ? "success" : "fail")   << endl;
    cout << "serializeDataTable3(): " << (serializeDataTable3() ? "success" : "fail")   << endl;
    cout << "serializeDataTable4(): " << (serializeDataTable4() ? "success" : "fail")   << endl;
    cout << "serializeDataTable5(): " << (serializeDataTable5() ? "success" : "fail")   << endl;
    cout << "serializeDataTable6(): " << (serializeDataTable6() ? "success" : "fail")   << endl;
    cout << "BSplines:                                                              "   << endl;
    cout << "serializeBSpline1():   " << (serializeBSpline1()   ? "success" : "fail")   << endl;
    cout << "serializePSpline1():   " << (serializePSpline1()   ? "success" : "fail")   << endl;
    cout << "serializeRBF1():       " << (serializeRBF1()       ? "success" : "fail")   << endl;

    return 0;
}