#ifndef TENSORINDEX_H
#define TENSORINDEX_H

// Class for transformation between tensor and vector indices
// Assuming zero-index
// Maybe rename to multi-index?
class TensorIndex {

public:

    // Constructor
    TensorIndex(std::vector<int> tensor_size)
        : tensor_size(tensor_size)
    {
        // Check tensor size and calculate required vector size
        vector_size = 1;
        assert(tensor_size.size() > 0);

        for (unsigned int i = 0; i < tensor_size.size(); i++)
        {
            assert(tensor_size.at(i) > 0);
            vector_size *= tensor_size.at(i);
        }
    }

    // Convert from N-D to 1-D
    int tensorToVectorIndex(const std::vector<int> ti) const
    {
        assert(checkTensorIndex(ti));

        // Calculate vector index vi
        int vi = ti.at(0);
        for (unsigned int i = 1; i < tensor_size.size(); i++)
        {
            vi = vi*tensor_size.at(i) + ti.at(i);
        }
        return vi;
    }

    // Convert from 1-D to N-D
    std::vector<int> vectorToTensorIndex(const int vi) const
    {
        int size = vector_size;
        int temp = vi;
        std::vector<int> ti;

        for (unsigned int i = 0; i < tensor_size.size(); i++)
        {
            size = size/tensor_size.at(i);
            ti.push_back(floor(temp / size));
            temp = temp - size*floor(temp / size);
        }

        //assert(tensorToVectorIndex(ti) == vi); // Check

        return ti;
    }

    // Check tensor indices
    bool checkTensorIndex(const std::vector<int> ti) const
    {
        for (unsigned int i = 0; i < ti.size(); i++)
        {
            if (ti.at(i) < 0 || ti.at(i) >= tensor_size.at(i))
            {
                return false;
            }
        }
        return true;
    }

    // Return vector size
    int vectorSize() const { return vector_size; }

    bool insideBounds(const std::vector<int> &t, const std::vector<int> &tlb, const std::vector<int> &tub ) const
    {
        for (unsigned int dim = 0; dim < t.size(); dim++)
        {
            if ( t.at(dim) < tlb.at(dim) || tub.at(dim) < t.at(dim) )
            {
                return false;
            }
        }
        return true;
    }// check  tlb <= t <= tub

    std::vector<int> tensorIndexSubtractBias(const std::vector<int> &t, const std::vector<int> &bias) const
    {
        std::vector<int> ret;
        for (unsigned int dim = 0; dim < t.size(); dim++)
        {
            ret.push_back( t.at(dim) - bias.at(dim));
        }
        return ret;
    } // subtract bias between two tensor representations

private:
    int vector_size;
    std::vector<int> tensor_size;
    TensorIndex();
    TensorIndex& operator = (TensorIndex const& assign);
};

#endif // TENSORINDEX_H
