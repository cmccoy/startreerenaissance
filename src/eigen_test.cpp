#include <atomic>
#include <vector>
#include <iostream>
#include <boost/math/tools/minima.hpp>

#include "gtr.hpp"

using namespace gtr;
using boost::math::tools::brent_find_minima;

double star_likelihood(const GTRModel& model,
                       const std::vector<Sequence>& sequences)
{
    double result = 0.0;

#pragma omp parallel for reduction(+:result)
    for(size_t i = 0; i < sequences.size(); i++) {
        result += model.logLikelihood(sequences[i]);
    }
    return result;
}

void estimate_branch_lengths(const GTRModel& model,
                             std::vector<Sequence>& sequences)
{
#pragma omp parallel for
    for(size_t i = 0; i < sequences.size(); i++) {
        Sequence& s = sequences[i];
        auto f = [&model, &s](const double d) {
            s.distance = d;
            return - model.logLikelihood(s);
        };
        boost::uintmax_t max_iter = 50;
        s.distance = brent_find_minima(f, 1e-7, 1.0, 8, max_iter).second;
    }
}

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
