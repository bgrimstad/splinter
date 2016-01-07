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
#include "datapoint.h"
#include <datatable.h>
#include <bspline.h>
#include <bsplinebasis.h>
#include <bsplinebasis1d.h>
#include "rbfnetwork.h"
#include <polynomial.h>

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
    //assert(ifs);

    stream.clear();
    // Convert from char to uint_8 vector
    for (const char& byte : result)
        stream.push_back((uint8_t)byte);

    read = stream.cbegin();
}

/*
 * get_size implementations
 */

size_t Serializer::get_size(const DenseVector &obj)
{
    size_t size = sizeof(obj.size());
    size += obj.size() * sizeof(double);
    return size;
}

size_t Serializer::get_size(const DataPoint &obj)
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

size_t Serializer::get_size(const RBFNetwork &obj)
{
    return get_size(obj.samples)
           + get_size(obj.normalized)
           + get_size(obj.precondition)
           + get_size(obj.numVariables)
           + get_size(obj.numSamples)
           + get_size(obj.type)
           + get_size(obj.coefficients);
}

size_t Serializer::get_size(const Polynomial &obj)
{
    return get_size(obj.numVariables)
            + get_size(obj.coefficients)
            + get_size(obj.degrees);
}

/*
 * _serialize implementations
 */

void Serializer::_serialize(const DenseVector &obj)
{
    _serialize(obj.size());
    for (size_t i = 0; i < obj.size(); ++i) {
        _serialize(obj(i));
    }
}

void Serializer::_serialize(const DataPoint &obj)
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

void Serializer::_serialize(const RBFNetwork &obj)
{
    _serialize(obj.samples);
    _serialize(obj.normalized);
    _serialize(obj.precondition);
    _serialize(obj.numVariables);
    _serialize(obj.numSamples);
    _serialize(obj.type);
    _serialize(obj.coefficients);
}

void Serializer::_serialize(const Polynomial &obj)
{
    _serialize(obj.numVariables);
    _serialize(obj.coefficients);
    _serialize(obj.degrees);
}

/*
 * deserialize implementations
 */

void Serializer::deserialize(DenseVector &obj)
{
    size_t length;
    deserialize(length);
    obj.resize(length);
    for (size_t i = 0; i < length; ++i) {
        deserialize(obj(i));
    }
}

void Serializer::deserialize(DataPoint &obj)
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

void Serializer::deserialize(RBFNetwork &obj)
{
    deserialize(obj.samples);
    deserialize(obj.normalized);
    deserialize(obj.precondition);
    deserialize(obj.numVariables);
    deserialize(obj.numSamples);
    deserialize(obj.type);
    if (obj.type == RBFType::THIN_PLATE_SPLINE)
    {
        obj.fn = std::shared_ptr<RBF>(new ThinPlateSpline());
    }
    else if (obj.type == RBFType::MULTIQUADRIC)
    {
        obj.fn = std::shared_ptr<RBF>(new Multiquadric());
    }
    else if (obj.type == RBFType::INVERSE_QUADRIC)
    {
        obj.fn = std::shared_ptr<RBF>(new InverseQuadric());
    }
    else if (obj.type == RBFType::INVERSE_MULTIQUADRIC)
    {
        obj.fn = std::shared_ptr<RBF>(new InverseMultiquadric());
    }
    else if (obj.type == RBFType::GAUSSIAN)
    {
        obj.fn = std::shared_ptr<RBF>(new Gaussian());
    }
    else
    {
        obj.fn = std::shared_ptr<RBF>(new ThinPlateSpline());
    }
    deserialize(obj.coefficients);
}

void Serializer::deserialize(Polynomial &obj)
{
    deserialize(obj.numVariables);
    deserialize(obj.coefficients);
    deserialize(obj.degrees);
}

} // namespace SPLINTER
