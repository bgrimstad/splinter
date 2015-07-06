/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_SERIALIZER_H
#define SPLINTER_SERIALIZER_H

#include <string>
#include <vector>
#include <iostream>
#include <generaldefinitions.h>
#include <set>
#include <stdint.h>

namespace SPLINTER
{

class DataSample;
class DataTable;
class BSpline;
class BSplineBasis;
class BSplineBasis1D;
class PSpline;
class RadialBasisFunction;
class PolynomialRegression;

class Serializer {
public:
    Serializer();
    Serializer(std::string fileName);

    // Serialize obj into the internal stream
    template <class T>
    void serialize(const T &obj);

    template <class T>
    void deserialize(T &obj);

    template <class T>
    void deserialize(std::vector<T> &obj);

    template <class T>
    void deserialize(std::set<T> &obj);

    template <class T>
    void deserialize(std::multiset<T> &obj);

    template <class T>
    void deserialize(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &obj);

    void deserialize(DataSample &obj);
    void deserialize(DataTable &obj);
    void deserialize(BSpline &obj);
    void deserialize(BSplineBasis &obj);
    void deserialize(BSplineBasis1D &obj);
    void deserialize(PSpline &obj);
    void deserialize(RadialBasisFunction &obj);
    void deserialize(PolynomialRegression &obj);

    // Save the serialized stream to fileName
    void saveToFile(std::string fileName);

    // Load fileName into the internal stream
    void loadFromFile(std::string fileName);

    virtual ~Serializer() {};

protected:
    template <class T>
    size_t get_size(const T &obj);

    template <class T>
    size_t get_size(const std::vector<T> &obj);

    template <class T>
    size_t get_size(const std::set<T> &obj);

    template <class T>
    size_t get_size(const std::multiset<T> &obj);

    template <class T>
    size_t get_size(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &obj);

    size_t get_size(const DataSample &obj);
    size_t get_size(const DataTable &obj);
    size_t get_size(const BSpline &obj);
    size_t get_size(const BSplineBasis &obj);
    size_t get_size(const BSplineBasis1D &obj);
    size_t get_size(const PSpline &obj);
    size_t get_size(const RadialBasisFunction &obj);
    size_t get_size(const PolynomialRegression &obj);

    template <class T>
    void _serialize(const T &obj);

    template <class T>
    void _serialize(const std::vector<T> &obj);

    template <class T>
    void _serialize(const std::set<T> &obj);

    template <class T>
    void _serialize(const std::multiset<T> &obj);

    template <class T>
    void _serialize(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &obj);

    void _serialize(const DataSample &obj);
    void _serialize(const DataTable &obj);
    void _serialize(const BSpline &obj);
    void _serialize(const BSplineBasis &obj);
    void _serialize(const BSplineBasis1D &obj);
    void _serialize(const PSpline &obj);
    void _serialize(const RadialBasisFunction &obj);
    void _serialize(const PolynomialRegression &obj);

private:
    typedef std::vector<uint8_t> StreamType;
    StreamType stream;

    // Where we are when serializing
    StreamType::iterator write;

    // Where we are when deserializing
    StreamType::const_iterator read;

};


template <class T>
void Serializer::serialize(const T &obj)
{
    // We can't set write to stream.end() here because the call
    // to stream.resize() below may invalidate iterators
    int writeIndex = stream.size();

    // Increase the size of the stream so it can hold the object
    stream.resize(stream.size() + get_size(obj));

    write = stream.begin() + writeIndex;

    _serialize(obj);
}

template <class T>
void Serializer::_serialize(const T &obj)
{
    // Get a uint8_t pointer to the object, so we can copy it into the stream
    auto objPtr = reinterpret_cast<const uint8_t *>(&obj);

    std::copy(objPtr, objPtr + sizeof(T), write);

    write += sizeof(T);
}

template <class T>
void Serializer::deserialize(T &obj)
{
    if(read + sizeof(T) > stream.cend()) {
        throw Exception("Serializer::deserialize: Stream is missing bytes!");
    }

    auto objPtr = reinterpret_cast<uint8_t *>(&obj);

    // Copy the data into val
    std::copy(read, read + sizeof(T), objPtr);

    read += sizeof(T);
}

template <class T>
size_t Serializer::get_size(const T &obj)
{
    return sizeof(T);
}

/*
 * get_size specializations
 */
template <class T>
size_t Serializer::get_size(const std::vector<T> &obj)
{
    size_t size = sizeof(size_t);
    for(auto &elem : obj) {
        size += get_size(elem);
    }

    return size;
}

template <class T>
size_t Serializer::get_size(const std::set<T> &obj)
{
    size_t size = sizeof(size_t);
    for(auto &elem : obj) {
        size += get_size(elem);
    }

    return size;
}

template <class T>
size_t Serializer::get_size(const std::multiset<T> &obj)
{
    size_t size = sizeof(size_t);
    for(auto &elem : obj) {
        size += get_size(elem);
    }

    return size;
}

template <class T>
size_t Serializer::get_size(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &obj)
{
    size_t size = sizeof(obj.rows());
    size += sizeof(obj.cols());
    size += obj.rows() * obj.cols() * sizeof(T);
    return size;
}


/*
 * _serialize specializations
 */
template <class T>
void Serializer::_serialize(const std::vector<T> &obj)
{
    _serialize(obj.size());
    for(auto &elem : obj)
    {
        _serialize(elem);
    }
}

template <class T>
void Serializer::_serialize(const std::set<T> &obj)
{
    _serialize(obj.size());
    for(auto &elem : obj)
    {
        _serialize(elem);
    }
}

template <class T>
void Serializer::_serialize(const std::multiset<T> &obj)
{
    _serialize(obj.size());
    for(auto &elem : obj)
    {
        _serialize(elem);
    }
}


template <class T>
void Serializer::_serialize(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &obj)
{
    // Store the number of matrix rows and columns first
    _serialize(obj.rows());
    _serialize(obj.cols());
    // Store the matrix elements
    for (size_t i = 0; i < obj.rows(); ++i) {
        for (size_t j = 0; j < obj.cols(); ++j) {
            _serialize(obj(i,j));
        }
    }
}


/*
 * deserialize specializations
 */
template <class T>
void Serializer::deserialize(std::vector<T> &obj)
{
    size_t size; deserialize(size);
    obj.resize(size);

    for(auto &elem : obj)
    {
        deserialize(elem);
    }
}

template <class T>
void Serializer::deserialize(std::set<T> &obj)
{
    size_t size; deserialize(size);

    T elem;
    for(size_t i = 0; i < size; ++i)
    {
        deserialize(elem);
        obj.insert(elem);
    }
}

template <class T>
void Serializer::deserialize(std::multiset<T> &obj)
{
    size_t size; deserialize(size);

    T elem;
    for(int i = 0; i < size; i++)
    {
        deserialize(elem);
        obj.insert(elem);
    }
}

template <class T>
void Serializer::deserialize(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &obj)
{
    // Retrieve the number of rows
    size_t rows; deserialize(rows);
    size_t cols; deserialize(cols);

    obj.resize(rows, cols);

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            deserialize(obj(i, j));
        }
    }
}

} // namespace SPLINTER

#endif // SPLINTER_SERIALIZER_H