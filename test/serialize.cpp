#include <datatable.h>
#include <datasample.h>
#include <serialize.h>
#include <iostream>

using namespace MultivariateSplines;

// Checks if a is within margin of b
bool equalsWithinRange(double a, double b, double margin = 0.0)
{
    return b - margin <= a && a <= b + margin;
}

bool is_identical(DataTable &a, DataTable &b)
{
    if (a.getNumVariables() != b.getNumVariables())
        return false;

    auto ait = a.cbegin(), bit = b.cbegin();
    for (; ait != a.cend() && bit != b.cend(); ait++, bit++)
    {
        for (unsigned int i = 0; i < a.getNumVariables(); i++)
        {
//            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getX().at(i) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getX().at(i) << " ";
            if (!equalsWithinRange(ait->getX().at(i), bit->getX().at(i)))
                return false;
        }

//            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getY().at(j) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getY().at(j) << " ";
        if (!equalsWithinRange(ait->getY(), bit->getY()))
            return false;
//        std::cout << std::endl;
    }

//    std::cout << "Finished comparing samples..." << std::endl;

    return ait == a.cend() && bit == b.cend();
}

bool test()
{
    DataTable table;

    auto x = std::vector<double>(4);
    double y;
    for (double i = -0.0001; i <= 0.0001; i += 0.000001)
    {
        for (double j = -0.01; j <= 0.01; j += 0.001)
        {
            for (double k = -0.01; k <= 0.01; k += 0.001)
            {
                for (double l = -100000.0; l < 0.0; l += 13720.0)
                {
                    x.at(0) = i;
                    x.at(1) = j;
                    x.at(2) = k;
                    x.at(3) = l;
                    y = i * j;
                    table.addSample(x, y);
                }
            }
        }
    }

    std::cout << "Size of serialized table: " << get_size(table) << std::endl;
    exit(1);

    StreamType stream;
    serialize(table, stream);
    save_to_file("test.datatable", stream);

    DataTable deserializedTable = deserialize<DataTable>(load_from_file("test.datatable"));

    return is_identical(table, deserializedTable);
}

int main() {
    if(test()) {
        std::cout << "Equal!" << std::endl;
    }

    return 0;
}
