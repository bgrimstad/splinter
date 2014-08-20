#include <vector>
#include <sorteddatatable.h>
#include <sstream>
#include <cmath> // abs


// Checks if a is within tolerance % of b
bool equalsWithin(double a, double b, double tolerance = 0.01)
{
    double bAbs = abs(b);
    double bMin = b - tolerance * bAbs;
    double bMax = b + tolerance * bAbs;
    return bMin <= a && a <= bMax;
}

// Checks if a is within margin of b
bool equalsWithinRange(double a, double b, double margin = 0.00001)
{
    return b - margin <= a && a <= b + margin;
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
//            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getX().at(i) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getX().at(i) << " ";
            if(!equalsWithinRange(ait->getX().at(i), bit->getX().at(i)))
                return false;
        }
        for(int j = 0; j < a.getDimY(); j++)
        {
//            std::cout << std::setprecision(SAVE_DOUBLE_PRECISION) << ait->getY().at(j) << " == " << std::setprecision(SAVE_DOUBLE_PRECISION) << bit->getY().at(j) << " ";
            if(!equalsWithinRange(ait->getY().at(j), bit->getY().at(j)))
                return false;
        }
//        std::cout << std::endl;
    }

//    std::cout << "Finished comparing samples..." << std::endl;

    return ait == a.cend() && bit == b.cend();
}


bool test1()
{
    SortedDataTable table;

    auto x = std::vector<double>(1);
    auto y = std::vector<double>(1);
    for(double i = -0.3; i <= 0.3; i += 0.04)
    {
        x.at(0) = i;
        y.at(0) = 2 * i;
        table.addSample(x, y);
    }

    table.saveDataTable("test1.datatable");

    SortedDataTable loadedTable;
    loadedTable.loadDataTable("test1.datatable");

    return is_identical(table, loadedTable);
}

bool test2()
{
    SortedDataTable table;

    auto x = std::vector<double>(2);
    auto y = std::vector<double>(1);
    for(double i = -0.3; i <= 0.3; i += 0.04)
    {
        for(double j = -0.4; j <= 1.0; j += 0.08)
        {
            x.at(0) = i;
            x.at(1) = j;
            y.at(0) = i * j;
            table.addSample(x, y);
        }
    }

    table.saveDataTable("test2.datatable");

    SortedDataTable loadedTable;
    loadedTable.loadDataTable("test2.datatable");

    return is_identical(table, loadedTable);
}

bool test3()
{
    SortedDataTable table;

    auto x = std::vector<double>(2);
    auto y = std::vector<double>(2);
    for(double i = -0.3; i <= 0.3; i += 0.04)
    {
        for(double j = -0.4; j <= 1.0; j += 0.03)
        {
            x.at(0) = i;
            x.at(1) = j;
            y.at(0) = i * j;
            y.at(1) = 2 * i*i - j;
            table.addSample(x, y);
        }
    }

    table.saveDataTable("test3.datatable");

    SortedDataTable loadedTable;
    loadedTable.loadDataTable("test3.datatable");

    return is_identical(table, loadedTable);
}

bool test4()
{
    SortedDataTable table;

    auto x = std::vector<double>(3);
    auto y = std::vector<double>(2);
    for(double i = -0.0001; i <= 0.0001; i += 0.000001)
    {
        for(double j = -0.01; j <= 0.01; j += 0.001)
        {
            for(double k = -0.01; k <= 0.01; k += 0.001)
            {
                x.at(0) = i;
                x.at(1) = j;
                x.at(2) = k;
                y.at(0) = i * j;
                y.at(1) = i * j * k;
                table.addSample(x, y);
            }
        }
    }

    table.saveDataTable("test4.datatable");

    SortedDataTable loadedTable;
    loadedTable.loadDataTable("test4.datatable");

    return is_identical(table, loadedTable);
}

bool test5()
{
    SortedDataTable table;

    auto x = std::vector<double>(4);
    auto y = std::vector<double>(2);
    for(double i = -0.0001; i <= 0.0001; i += 0.000001)
    {
        for(double j = -0.01; j <= 0.01; j += 0.001)
        {
            for(double k = -0.01; k <= 0.01; k += 0.001)
            {
                for(double l = -100000.0; l < 0.0; l += 13720.0)
                {
                    x.at(0) = i;
                    x.at(1) = j;
                    x.at(2) = k;
                    x.at(3) = l;
                    y.at(0) = i * j;
                    y.at(1) = i * j * k;
                    y.at(1) = i * j * k * l;
                    table.addSample(x, y);
                }
            }
        }
    }

    table.saveDataTable("test5.datatable");

    SortedDataTable loadedTable;
    loadedTable.loadDataTable("test5.datatable");

    return is_identical(table, loadedTable);
}


int main(int argc, char **argv)
{
    std::cout << "test1(): " << (test1() ? "success" : "fail") << std::endl;
    std::cout << "test2(): " << (test2() ? "success" : "fail") << std::endl;
    std::cout << "test3(): " << (test3() ? "success" : "fail") << std::endl;
    std::cout << "test4(): " << (test4() ? "success" : "fail") << std::endl;
    std::cout << "test5(): " << (test5() ? "success" : "fail") << std::endl;
}
