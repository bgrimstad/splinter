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
#include <utils/test_function_collection.h>
#include <utils/timer.h>


int main(int argc, char *argv[])
{
    setup_test_functions();

    SPLINTER::Timer tim;
    tim.start();

    int result = Catch::Session().run(argc, argv);

    tim.stop();
    std::cout << "Test run time: " << tim.get_milli_seconds() << " (ms)" << std::endl;

    tearDownTestFunctions();

    return result;
}
