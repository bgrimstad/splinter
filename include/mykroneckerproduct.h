/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef MS_MYKRONECKERPRODUCT_H
#define MS_MYKRONECKERPRODUCT_H

#include "generaldefinitions.h"

namespace MultivariateSplines
{

void myKroneckerProduct(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &AB);

} // namespace MultivariateSplines

#endif // MS_MYKRONECKERPRODUCT_H
