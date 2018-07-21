/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <utils/op_overloads.h>
#include <Catch.h>
#include <utilities.h>
#include <utils/test_utils.h>
#include <limits>
#include <iostream>

using namespace SPLINTER;


#define COMMON_TAGS "[general][utilities]"


TEST_CASE("Eigen types to standard C++ conversion", COMMON_TAGS)
{
    std::vector<double> originalVec =
    {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            0,-1,-2,-3,-4,-5,-6,-7,-8,-9,
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::min()
    };
    std::vector<std::vector<double>> originalVecVec =
    {
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
            {0,-1,-2,-3,-4,-5,-6,-7,-8,-9},
            {std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()},
            {-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity(),-std::numeric_limits<double>::infinity()},
            {std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max()},
            {std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min(),std::numeric_limits<double>::min()}
    };

    bool testVectorOk = true;
    size_t cols = originalVecVec.at(0).size();
    for (size_t i = 0; i < originalVecVec.size(); ++i)
    {
        testVectorOk = testVectorOk && (originalVecVec.at(i).size() == cols);
    }

    {
        INFO("Test vector is not rectangular!");
        REQUIRE(testVectorOk);
    }

    std::vector<double> emptyVec;
    std::vector<std::vector<double>> emptyVecVec;
    DenseVector emptyDenseVec;
    DenseMatrix emptyDenseMat;


    SECTION("Converting between standard vector and DenseVector works correctly")
    {
        CHECK(originalVec == std_to_eig_vec(originalVec));
        CHECK(originalVec == eig_to_std_vec(std_to_eig_vec(originalVec)));
        CHECK(emptyVec == eig_to_std_vec(emptyDenseVec));
        CHECK(emptyVec == eig_to_std_vec(std_to_eig_vec(emptyVec)));
    }

    SECTION("Converting between standard vector<vector<double>> and DenseMatrix works correctly")
    {
        CHECK(originalVecVec == std_to_eig_mat(originalVecVec));
        CHECK(originalVecVec == eig_to_std_mat(std_to_eig_mat(originalVecVec)));
        CHECK(emptyVecVec == std_to_eig_mat(emptyVecVec));
        CHECK(emptyVecVec == eig_to_std_mat(std_to_eig_mat(emptyVecVec)));
    }
}
