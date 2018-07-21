/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_TESTFUNCTIONS_H
#define SPLINTER_TESTFUNCTIONS_H

#include <utils/test_function.h>

extern std::vector<std::vector<SPLINTER::TestFunction *>> testFunctions;

void setup_test_functions();

void tearDownTestFunctions();

#endif // SPLINTER_TESTFUNCTIONS_H
