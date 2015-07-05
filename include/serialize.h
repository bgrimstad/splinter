/**
* Copyright (c) 2013, Simone Pellegrini All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* - Redistributions of source code must retain the above copyright notice,
* this list of conditions and the following disclaimer.
*
* - Redistributions in binary form must reproduce the above copyright notice,
* this list of conditions and the following disclaimer in the documentation
* and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
* Modifications to original code by Bjarne Grimstad:
* - Removed support for tuples and strings (not needed here)
* - Added support for dense Eigen matrices
* - Added support for saving to and loading from files
*/

#ifndef SPLINTER_SERIALIZE_H
#define SPLINTER_SERIALIZE_H

#include <generaldefinitions.h>
#include <vector>
#include <Eigen/Dense>
#include <cassert>
#include <fstream>
#include <numeric>
#include <set>
#include <bspline.h>
#include <datasample.h>
#include <stdint.h>
#include <bsplinebasis.h>
#include <pspline.h>
#include <radialbasisfunction.h>

namespace SPLINTER
{

template <class T>
size_t get_size(const T& obj);

namespace detail
{
    template <class T>
    struct get_size_helper;

    /**
     * Specialization for Eigen's dense matrix. A matrix is represented
     * in the stream by storing the number of rows (mat.rows()),
     * the number of columns (mat.cols()), and the matrix elements by
     * iteration over rows, then columns.
     */
    template <class T>
    struct get_size_helper<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> {
        static size_t value(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& obj) {
            size_t size = sizeof(obj.rows());
            size += sizeof(obj.cols());
            size += obj.rows()*obj.cols()*sizeof(T);
            return size;
        }
    };

    /**
     * Specialization for DataSample.
     */
    template <>
    struct get_size_helper<DataSample> {
        static size_t value(const DataSample& sample) {
            return get_size(sample.getX()) + get_size(sample.getY());
        }
    };

    template <>
    struct get_size_helper<DataTable> {
        static size_t value(const DataTable& table) {
            return get_size(table.getAllowDuplicates())
                   + get_size(table.getAllowIncompleteGrid())
                   + get_size(table.getNumDuplicates())
                   + get_size(table.getNumVariables())
                   + get_size(table.getSamples())
                   + get_size(table.getGrid());
        }
    };

    template <>
    struct get_size_helper<BSpline> {
        static size_t value(const BSpline &bspline) {
            return get_size(bspline.getBasis())
                   + get_size(bspline.getKnotaverages())
                   + get_size(bspline.getCoefficients())
                   + get_size(bspline.getNumVariables());
        }
    };

    template <>
    struct get_size_helper<PSpline> {
        static size_t value(const PSpline &pspline) {
            return get_size(pspline.getLambda()) + get_size(*((BSpline *) &pspline));
        }
    };

    template <>
    struct get_size_helper<RadialBasisFunction> {
        static size_t value(const RadialBasisFunction &rbf) {
            return get_size(rbf.samples)
                   + get_size(rbf.normalized)
                   + get_size(rbf.precondition)
                   + get_size(rbf.dim)
                   + get_size(rbf.numSamples)
                   + get_size(rbf.type)
                   + get_size(rbf.weights);
        }
    };

    template <>
    struct get_size_helper<BSplineBasis> {
        static size_t value(const BSplineBasis &basis) {
            return get_size(basis.getBases()) + get_size(basis.getNumVariables());
        }
    };

    template <>
    struct get_size_helper<BSplineBasis1D> {
        static size_t value(const BSplineBasis1D &basis) {
            return get_size(basis.getDegree()) + get_size(basis.getKnots()) + get_size(basis.getTargetNumBasisFunctions());
        }
    };

    /**
    * Specialization for generic std::vector<T>. A vector is represented
    * in the stream by storing the number of elements (v.size()) followed
    * by the serialized elements (notice that we can have nested structures).
    */
    template <class T>
    struct get_size_helper<std::vector<T>> {
        static size_t value(const std::vector<T>& obj) {
            return std::accumulate(obj.begin(), obj.end(), sizeof(size_t),
                [](const size_t& acc, const T& cur) { return acc+get_size(cur); });
        }
    };

    /**
    * Specialization for generic std::multiset<T>. A multiset is represented
    * in the stream by storing the number of elements (m.size()) followed
    * by the serialized elements (notice that we can have nested structures).
    */
    template <class T>
    struct get_size_helper<std::multiset<T>> {
        static size_t value(const std::multiset<T>& obj) {
            return std::accumulate(obj.begin(), obj.end(), sizeof(size_t),
                [](const size_t& acc, const T& cur) { return acc+get_size(cur); });
        }
    };

    /**
    * Specialization for generic std::set<T>. A multiset is represented
    * in the stream by storing the number of elements (s.size()) followed
    * by the serialized elements (notice that we can have nested structures).
    */
    template <class T>
    struct get_size_helper<std::set<T>> {
        static size_t value(const std::set<T>& obj) {
            return std::accumulate(obj.begin(), obj.end(), sizeof(size_t),
                [](const size_t& acc, const T& cur) { return acc+get_size(cur); });
        }
    };

    /**
    * Specialization for any remaining type, simply use the value of
    * sizeof()
    */
    template <class T>
    struct get_size_helper {
        static size_t value(const T& obj) { return sizeof(T); }
    };
} // end detail namespace

template <class T>
inline size_t get_size(const T& obj) {
    return detail::get_size_helper<T>::value(obj);
}

namespace detail
{
    template <class T>
    struct serialize_helper;

    template <class T>
    void serializer(const T& obj, StreamType::iterator&);

    template <class T>
    struct serialize_helper<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> {
        static void apply(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& obj, StreamType::iterator& res) {
            // Store the number of matrix rows and columns first
            serializer(obj.rows(), res);
            serializer(obj.cols(), res);
            // Store the matrix elements
            for (int i = 0; i < obj.rows(); ++i) {
                for (int j = 0; j < obj.cols(); ++j) {
                    serializer(obj(i,j), res);
                }
            }
        }
    };

    /**
     * Serialization of DataSample.
     */
    template <>
    struct serialize_helper<DataSample> {
        static void apply(const DataSample& sample, StreamType::iterator& res) {
            serializer(sample.getX(), res);
            serializer(sample.getY(), res);
        }
    };

    template <>
    struct serialize_helper<DataTable> {
        static void apply(const DataTable &table, StreamType::iterator &res) {
            serializer(table.getAllowDuplicates(), res);
            serializer(table.getAllowIncompleteGrid(), res);
            serializer(table.getNumDuplicates(), res);
            serializer(table.getNumVariables(), res);
            serializer(table.getSamples(), res);
            serializer(table.getGrid(), res);
        }
    };

    template <class T>
    struct serialize_helper<std::vector<T>> {
        static void apply(const std::vector<T>& obj, StreamType::iterator& res) {
            // Store the number of vector elements first
            serializer(obj.size(), res);
            for (const auto& cur : obj) { serializer(cur, res); }
        }
    };

    template <class T>
    struct serialize_helper<std::multiset<T>> {
        static void apply(const std::multiset<T>& obj, StreamType::iterator& res) {
            // Store the number of multiset elements first
            serializer(obj.size(), res);
            for (const auto& cur : obj) { serializer(cur, res); }
        }
    };

    template <class T>
    struct serialize_helper<std::set<T>> {
        static void apply(const std::set<T>& obj, StreamType::iterator& res) {
            // Store the number of set elements first
            serializer(obj.size(), res);
            for (const auto& cur : obj) { serializer(cur, res); }
        }
    };

    template <>
    struct serialize_helper<BSpline> {
        static void apply(const BSpline &obj, StreamType::iterator &res)
        {
            serializer(obj.getBasis(), res);

            serializer(obj.getKnotaverages(), res);

            serializer(obj.getCoefficients(), res);

            serializer(obj.getNumVariables(), res);
        }
    };

    template <>
    struct serialize_helper<PSpline> {
        static void apply(const PSpline &obj, StreamType::iterator &res)
        {
            serializer(*(BSpline *) &obj, res);

            serializer(obj.getLambda(), res);
        }
    };

    template <>
    struct serialize_helper<BSplineBasis> {
        static void apply(const BSplineBasis &obj, StreamType::iterator &res)
        {
            // Store the bases
            serializer(obj.getBases(), res);

            serializer(obj.getNumVariables(), res);
        }
    };

    template <>
    struct serialize_helper<BSplineBasis1D> {
        static void apply(const BSplineBasis1D &obj, StreamType::iterator &res)
        {
            serializer(obj.getDegree(), res);

            serializer(obj.getKnots(), res);

            serializer(obj.getTargetNumBasisFunctions(), res);
        }
    };

    template <>
    struct serialize_helper<RadialBasisFunction> {
        static void apply(const RadialBasisFunction &obj, StreamType::iterator &res)
        {
            serializer(obj.samples, res);
            serializer(obj.normalized, res);
            serializer(obj.precondition, res);
            serializer(obj.dim, res);
            serializer(obj.numSamples, res);
            serializer(obj.type, res);
            serializer(obj.weights, res);
        }
    };

    template <class T>
    struct serialize_helper {
        static void apply(const T& obj, StreamType::iterator& res) {
            const uint8_t* ptr = reinterpret_cast<const uint8_t*>(&obj);
            // Copy the contents between ptr and ptr+sizeof(T) to res
            // The returned value is the element beyond the last element copied
            res = std::copy(ptr, ptr+sizeof(T), res);
        }
    };

    template <class T>
    inline void serializer(const T& obj, StreamType::iterator& res) {
        serialize_helper<T>::apply(obj, res);
    }
} // end detail namespace

template <class T>
inline void serialize(const T& obj, StreamType& res) {
    size_t offset = res.size();
    size_t size = get_size(obj);
    res.resize(res.size() + size);

    StreamType::iterator it = res.begin() + offset;
    detail::serializer(obj,it);
    assert(res.begin() + offset + size == it);
}

namespace detail
{
    template <class T>
    struct deserialize_helper;

    /**
    * Deserialization for integral types and POD data types.
    *
    * This is done by simply relying on the sizeof() operator of the object.
    * It is important that the datatype we are deserializing has a default
    * contructor, otherwise you need to provide a specialization for that type
    */
    template <class T>
    struct deserialize_helper {
        static T apply(StreamType::const_iterator& begin,
                       StreamType::const_iterator end) {
            assert(begin+sizeof(T)<=end && "Error: not enough bytes in the stream to deserialize type");
            T val;
            uint8_t* ptr = reinterpret_cast<uint8_t*>(&val);
            std::copy(begin, begin+sizeof(T), ptr);
            begin += sizeof(T);
            return val;
        }
    };

    /**
    * Deserialization for vector types.
    *
    * We read the first size_t value which contains the number of elements
    * and then recursively read each of them.
    */
    template <class T>
    struct deserialize_helper<std::vector<T>> {
        static std::vector<T> apply(StreamType::const_iterator& begin,
                                    StreamType::const_iterator end)
        {
            // Retrieve the number of elements
            size_t size = deserialize_helper<size_t>::apply(begin, end);

            std::vector<T> vect(size);
            for (size_t i = 0; i < size; ++i)
            {
                /**
                * Call the move-copy constructor so that the additional copy
                * is avoided
                */
                vect[i] = std::move(deserialize_helper<T>::apply(begin, end));
            }
            return vect;
        }
    };

    /**
    * Deserialization for multiset types.
    *
    * We read the first size_t value which contains the number of elements
    * and then recursively read each of them.
    */
    template <class T>
    struct deserialize_helper<std::multiset<T>> {
        static std::multiset<T> apply(StreamType::const_iterator& begin,
                                    StreamType::const_iterator end)
        {
            // Retrieve the number of elements
            size_t size = deserialize_helper<size_t>::apply(begin, end);

            std::multiset<T> multiset;
            for (size_t i = 0; i < size; ++i)
            {
                multiset.insert(deserialize_helper<T>::apply(begin, end));
            }

            return multiset;
        }
    };

    /**
    * Deserialization for set types.
    *
    * We read the first size_t value which contains the number of elements
    * and then recursively read each of them.
    */
    template <class T>
    struct deserialize_helper<std::set<T>> {
        static std::set<T> apply(StreamType::const_iterator& begin,
                                    StreamType::const_iterator end)
        {
            // Retrieve the number of elements
            size_t size = deserialize_helper<size_t>::apply(begin, end);

            std::set<T> set;
            for (size_t i = 0; i < size; ++i)
            {
                set.insert(deserialize_helper<T>::apply(begin, end));
            }
            return set;
        }
    };

    /**
    * Deserialization for Eigen dense matrices.
    *
    * We read the first size_t values for number of rows and columns,
    * then we read the matrix elements doubleby iterating over rows, then columns.
    */
    template <class T>
    struct deserialize_helper<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> {
        static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> apply(StreamType::const_iterator& begin,
                                    StreamType::const_iterator end)
        {
            // Retrieve the number of rows
            size_t rows = deserialize_helper<size_t>::apply(begin, end);
            size_t cols = deserialize_helper<size_t>::apply(begin, end);

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat(rows, cols);
            for (size_t i = 0; i < rows; ++i)
            {
                for (size_t j = 0; j < cols; ++j)
                {
                    /**
                    * Call the move-copy constructor so that the additional copy
                    * is avoided
                    */
                    mat(i,j) = std::move(deserialize_helper<T>::apply(begin, end));
                }

            }
            return mat;
        }
    };

    /**
    * Deserialization for DataSample.
    */
    template <>
    struct deserialize_helper<DataSample> {
        static DataSample apply(StreamType::const_iterator& begin, StreamType::const_iterator end) {
            auto x = deserialize_helper<std::vector<double>>::apply(begin, end);
            auto y = deserialize_helper<double>::apply(begin, end);

            DataSample sample(x, y);

            return sample;
        }
    };

    template <>
    struct deserialize_helper<DataTable> {
        static DataTable apply(StreamType::const_iterator& begin, StreamType::const_iterator end) {
            DataTable table;
            table._deserialize(begin, end);
            return table;
        }
    };

    template <>
    struct deserialize_helper<BSplineBasis1D> {
        static BSplineBasis1D apply(StreamType::const_iterator& begin, StreamType::const_iterator end)
        {
            auto degree = deserialize_helper<unsigned int>::apply(begin, end);

            auto knots = deserialize_helper<std::vector<double>>::apply(begin, end);

            auto targetNumBasisFunctions = deserialize_helper<unsigned int>::apply(begin, end);

            BSplineBasis1D basis1d(degree, knots, targetNumBasisFunctions);

            return basis1d;
        }
    };

    template <>
    struct deserialize_helper<BSplineBasis> {
        static BSplineBasis apply(StreamType::const_iterator& begin, StreamType::const_iterator end)
        {
            auto bases = deserialize_helper<std::vector<BSplineBasis1D>>::apply(begin, end);

            auto numVariables = deserialize_helper<unsigned int>::apply(begin, end);

            BSplineBasis basis(bases, numVariables);

            return basis;
        }
    };

    template <>
    struct deserialize_helper<BSpline> {
        static BSpline apply(StreamType::const_iterator& begin, StreamType::const_iterator end)
        {
            BSpline bspline;
            bspline._deserialize(begin, end);
            return bspline;
        }
    };

    template <>
    struct deserialize_helper<PSpline> {
        static PSpline apply(StreamType::const_iterator& begin, StreamType::const_iterator end)
        {
            PSpline pspline;
            pspline._deserialize(begin, end);
            return pspline;
        }
    };

} // end namespace detail

template <class T>
inline T deserialize(StreamType::const_iterator& begin,
                     const StreamType::const_iterator& end) {
    return detail::deserialize_helper<T>::apply(begin, end);
}

template <class T>
inline T deserialize(const StreamType& res) {
    StreamType::const_iterator it = res.cbegin();
    return deserialize<T>(it, res.cend());
}

/**
 * Save serialized stream to file
 */
inline void save_to_file(std::string filename, StreamType data)
{
    std::fstream fs(filename, std::fstream::out | std::fstream::binary);

    for (const auto& byte : data)
        fs << byte;
}

/**
 * Loads serialized stream from file
 */
inline StreamType load_from_file(std::string filename)
{
    // Open the file in binary mode at the end
    std::ifstream ifs(filename, std::ios::binary | std::ios::ate);

    if(!ifs.is_open()) {
        std::string error_message("serialize.h: load_from_file: Unable to open file \"");
        error_message.append(filename);
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

    // Convert from char to uint_8 vector
    StreamType sresult;
    for (const char& byte : result)
        sresult.push_back((uint8_t)byte);

    return sresult;
}

} // namespace SPLINTER

#endif // SPLINTER_SERIALIZE_H
