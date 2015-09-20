/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "serializer.h"
#include <fstream>
#include "definitions.h"
#include <datasample.h>
#include <datatable.h>
#include <bspline.h>
#include <bsplinebasis.h>
#include <bsplinebasis1d.h>
#include <bsplineapproximant.h>
#include <psplineapproximant.h>
#include <rbfapproximant.h>
#include <polynomial.h>
#include <polynomialapproximant.h>

namespace SPLINTER
{

Serializer::Serializer()
{
    stream = StreamType(0);
}

Serializer::Serializer(std::string fileName)
{
    stream = StreamType(0);
    loadFromFile(fileName);
}

void Serializer::saveToFile(std::string fileName)
{
    std::fstream fs(fileName, std::fstream::out | std::fstream::binary);

    for (const auto& byte : stream)
        fs << byte;
}

void Serializer::loadFromFile(std::string fileName)
{
    // Open the file in binary mode at the end
    std::ifstream ifs(fileName, std::ios::binary | std::ios::ate);

    if (!ifs.is_open()) {
        std::string error_message("Serializer::loadFromFile: Unable to open file \"");
        error_message.append(fileName);
        error_message.append("\" for deserializing.");
        throw Exception(error_message);
    }

    // Because we opened the file at the end, tellg() will give us the size of the file
    std::ifstream::pos_type pos = ifs.tellg();

    std::vector<char> result(pos);

    ifs.seekg(0, std::ios::beg);

    // http://www.cplusplus.com/reference/vector/vector/data/
    // Elements of the vector are guaranteed to be stored in a contiguous array,
    // *result.data() can therefore be treated as an array of the same type as the vector
    ifs.read(result.data(), pos);
    assert(ifs);

    stream.clear();
    // Convert from char to uint_8 vector
    for (const char& byte : result)
        stream.push_back((uint8_t)byte);

    read = stream.cbegin();
}

/*
 * get_size implementations
 */

size_t Serializer::get_size(const DataSample &obj)
{
    return get_size(obj.x) + get_size(obj.y);
}

size_t Serializer::get_size(const DataTable &obj)
{
    return get_size(obj.allowDuplicates)
           + get_size(obj.allowIncompleteGrid)
           + get_size(obj.numDuplicates)
           + get_size(obj.numVariables)
           + get_size(obj.samples)
           + get_size(obj.grid);
}

size_t Serializer::get_size(const BSpline &obj)
{
    return get_size(obj.basis)
           + get_size(obj.knotaverages)
           + get_size(obj.coefficients)
           + get_size(obj.numVariables);
}

size_t Serializer::get_size(const BSplineBasis &obj)
{
    return get_size(obj.bases)
           + get_size(obj.numVariables);
}

size_t Serializer::get_size(const BSplineBasis1D &obj)
{
    return get_size(obj.degree)
           + get_size(obj.knots)
           + get_size(obj.targetNumBasisfunctions);
}

size_t Serializer::get_size(const BSplineApproximant &obj)
{
    return get_size(obj.numVariables)
           + get_size(obj.bspline);
}

size_t Serializer::get_size(const PSplineApproximant &obj)
{
    return get_size(obj.numVariables)
           + get_size(obj.bspline)
           + get_size(obj.lambda);
}

size_t Serializer::get_size(const RBFApproximant &obj)
{
    return get_size(obj.samples)
           + get_size(obj.normalized)
           + get_size(obj.precondition)
           + get_size(obj.numVariables)
           + get_size(obj.numSamples)
           + get_size(obj.type)
           + get_size(obj.weights);
}

size_t Serializer::get_size(const Polynomial &obj)
{
    return get_size(obj.numVariables)
            + get_size(obj.coefficients)
            + get_size(obj.degrees);
}

size_t Serializer::get_size(const PolynomialApproximant &obj)
{
    return get_size(obj.numVariables)
            + get_size(obj.poly);
}

/*
 * _serialize implementations
 */

void Serializer::_serialize(const DataSample &obj)
{
    _serialize(obj.x);
    _serialize(obj.y);
}

void Serializer::_serialize(const DataTable &obj)
{
    _serialize(obj.allowDuplicates);
    _serialize(obj.allowIncompleteGrid);
    _serialize(obj.numDuplicates);
    _serialize(obj.numVariables);
    _serialize(obj.samples);
    _serialize(obj.grid);
}

void Serializer::_serialize(const BSpline &obj)
{
    _serialize(obj.basis);
    _serialize(obj.knotaverages);
    _serialize(obj.coefficients);
    _serialize(obj.numVariables);
}

void Serializer::_serialize(const BSplineBasis &obj)
{
    _serialize(obj.bases);
    _serialize(obj.numVariables);
}

void Serializer::_serialize(const BSplineBasis1D &obj)
{
    _serialize(obj.degree);
    _serialize(obj.knots);
    _serialize(obj.targetNumBasisfunctions);
}

void Serializer::_serialize(const BSplineApproximant &obj)
{
    _serialize(obj.numVariables);
    _serialize(obj.bspline);
}

void Serializer::_serialize(const PSplineApproximant &obj)
{
    _serialize(obj.numVariables);
    _serialize(obj.bspline);
    _serialize(obj.lambda);
}

void Serializer::_serialize(const RBFApproximant &obj)
{
    _serialize(obj.samples);
    _serialize(obj.normalized);
    _serialize(obj.precondition);
    _serialize(obj.numVariables);
    _serialize(obj.numSamples);
    _serialize(obj.type);
    _serialize(obj.weights);
}

void Serializer::_serialize(const Polynomial &obj)
{
    _serialize(obj.numVariables);
    _serialize(obj.coefficients);
    _serialize(obj.degrees);
}

void Serializer::_serialize(const PolynomialApproximant &obj)
{
    _serialize(obj.numVariables);
    _serialize(obj.poly);
}

/*
 * deserialize implementations
 */

void Serializer::deserialize(DataSample &obj)
{
    deserialize(obj.x);
    deserialize(obj.y);
}

void Serializer::deserialize(DataTable &obj)
{
    deserialize(obj.allowDuplicates);
    deserialize(obj.allowIncompleteGrid);
    deserialize(obj.numDuplicates);
    deserialize(obj.numVariables);
    deserialize(obj.samples);
    deserialize(obj.grid);
}

void Serializer::deserialize(BSpline &obj)
{
    deserialize(obj.basis);
    deserialize(obj.knotaverages);
    deserialize(obj.coefficients);
    deserialize(obj.numVariables);
}

void Serializer::deserialize(BSplineBasis &obj)
{
    deserialize(obj.bases);
    deserialize(obj.numVariables);
}

void Serializer::deserialize(BSplineBasis1D &obj)
{
    deserialize(obj.degree);
    deserialize(obj.knots);
    deserialize(obj.targetNumBasisfunctions);
}

void Serializer::deserialize(BSplineApproximant &obj)
{
    deserialize(obj.numVariables);
    deserialize(obj.bspline);
}

void Serializer::deserialize(PSplineApproximant &obj)
{
    deserialize(obj.numVariables);
    deserialize(obj.bspline);
    deserialize(obj.lambda);
}

void Serializer::deserialize(RBFApproximant &obj)
{
    deserialize(obj.samples);
    deserialize(obj.normalized);
    deserialize(obj.precondition);
    deserialize(obj.numVariables);
    deserialize(obj.numSamples);
    deserialize(obj.type);
    if (obj.type == RBFType::THIN_PLATE_SPLINE)
    {
        obj.fn = std::shared_ptr<RBFTerm>(new ThinPlateSpline());
    }
    else if (obj.type == RBFType::MULTIQUADRIC)
    {
        obj.fn = std::shared_ptr<RBFTerm>(new Multiquadric());
    }
    else if (obj.type == RBFType::INVERSE_QUADRIC)
    {
        obj.fn = std::shared_ptr<RBFTerm>(new InverseQuadric());
    }
    else if (obj.type == RBFType::INVERSE_MULTIQUADRIC)
    {
        obj.fn = std::shared_ptr<RBFTerm>(new InverseMultiquadric());
    }
    else if (obj.type == RBFType::GAUSSIAN)
    {
        obj.fn = std::shared_ptr<RBFTerm>(new Gaussian());
    }
    else
    {
        obj.fn = std::shared_ptr<RBFTerm>(new ThinPlateSpline());
    }
    deserialize(obj.weights);
}

void Serializer::deserialize(Polynomial &obj)
{
    deserialize(obj.numVariables);
    deserialize(obj.coefficients);
    deserialize(obj.degrees);
}

void Serializer::deserialize(PolynomialApproximant &obj)
{
    deserialize(obj.numVariables);
    deserialize(obj.poly);
}

} // namespace SPLINTER
