#include <vector>
#include <iostream>

#include "gtr.hpp"

using namespace gtr;

int main()
{
    GTRParameters m;
    GTRModel model = m.buildModel();

    std::cout << "eigenvalues:\n";
    std::cout << model.decomp.eigenvalues() << '\n';

    std::cout << m.params << '\n';
    std::cout << m.pi << '\n';
    std::cout << m.buildQMatrix() << '\n';
    std::cout << "Distance = .1\n" << model.buildPMatrix(0.1) << '\n';
    return 0;
}
