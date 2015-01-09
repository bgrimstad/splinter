#include <datatable.h>
#include <datasample.h>
#include <serialize.h>
#include <iostream>

using namespace MultivariateSplines;

int main() {
    DataTable t;
    t.addSample(2.0, 1.0);

    std::cout << get_size(t) << std::endl;

    auto sample = *t.getSamples().begin();

    std::cout << get_size(sample) << std::endl;

    StreamType stream;
    serialize(t, stream);

    auto it = stream.cbegin();

    DataTable deserialized = deserialize<DataTable>(it, stream.cend());

    t.printSamples(std::cout);
    deserialized.printSamples(std::cout);

    return 0;
}
