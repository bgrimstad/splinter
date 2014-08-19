#include <vector>
#include <sorteddatatable.h>
#include <sstream>
#include <cmath> // abs


bool equalsWithin(double a, double b, double tolerance = 0.02)
{
    double bAbs = abs(b);
    double bMin = b - (tolerance / 2.0) * bAbs;
    double bMax = b + (tolerance / 2.0) * bAbs;
    return bMin <= a && a <= bMax;
}

bool equalsWithinRange(double a, double b, double range = 0.001)
{
    double bMin = b - (range / 2.0);
    double bMax = b + (range / 2.0);
    return bMin <= a && a <= bMax;
}


bool is_identical(SortedDataTable &a, SortedDataTable &b)
{
    if(a.getDimX() != b.getDimX() || a.getDimY() != b.getDimY())
        return false;

    auto ait = a.cbegin(), bit = b.cbegin();
    for(; ait != a.cend() && bit != b.cend(); ait++, bit++)
    {
        for(int i = 0; i < a.getDimX(); i++)
        {
            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getX().at(i) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getX().at(i) << " ";
            if(!equalsWithinRange(ait->getX().at(i), bit->getX().at(i)))
                return false;
        }
        for(int j = 0; j < a.getDimY(); j++)
        {
            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getY().at(j) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getY().at(j) << " ";
            if(!equalsWithinRange(ait->getY().at(j), bit->getY().at(j)))
                return false;
        }
        std::cout << std::endl;
    }

    std::cout << "Finished comparing samples.." << std::endl;

    return ait == a.cend() && bit == b.cend();
}


int main(int argc, char **argv)
{
    SortedDataTable table;

    auto x = std::vector<double>(2);
    auto y = std::vector<double>(1);
    for(double i = 0.0; i < 0.3; i += 0.1)
    {
        for(double j = 0.0; j < 0.3; j += 0.1)
        {
            x.at(0) = i;
            x.at(1) = j;
            y.at(0) = i * j;
            table.addSample(x, y);
        }
    }

    table.saveDataTable("test.datatable");

    SortedDataTable loadedTable;
    loadedTable.loadDataTable("test.datatable");

    if(is_identical(table, loadedTable))
        std::cout << "Success!" << std::endl;
    else
        std::cout << "Failure!" << std::endl;

    for(double i = 9.0; i <= 11.0; i += 0.01)
    {
        bool eq = equalsWithin(i, 10.0);
        if(eq)
            std::cout << i << std::endl;
    }
}
