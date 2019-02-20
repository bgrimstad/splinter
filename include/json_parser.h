/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_JSON_PARSER_H
#define SPLINTER_JSON_PARSER_H

#include <string>

namespace SPLINTER
{

class BSpline;
class DataTable;

void bspline_to_json(const BSpline &bspline, const std::string &filename);
BSpline bspline_from_json(const std::string &filename);

void datatable_to_json(const DataTable &data, const std::string &filename);
DataTable datatable_from_json(const std::string &filename);

} // namespace SPLINTER

#endif //SPLINTER_JSON_PARSER_H
