#ifndef STAROPTIM_STAROPTIM_H
#define STAROPTIM_STAROPTIM_H

#include <cstdlib>
#include <vector>

struct Sequence;

namespace bpp
{
class SubstitutionModel;
class DiscreteDistribution;
}

namespace star_optim
{

/// \brief calculate the star-tree likelihood
double starLikelihood(const bpp::SubstitutionModel&,
                      const bpp::DiscreteDistribution&,
                      const std::vector<Sequence>&);

/// \brief estimate branch lengths
void estimateBranchLengths(const bpp::SubstitutionModel&,
                           const bpp::DiscreteDistribution&,
                           std::vector<Sequence>&);

/// \brief Optimize the model & branch lengths distribution for a collection of sequences
size_t optimize(bpp::SubstitutionModel&,
                bpp::DiscreteDistribution&,
                std::vector<Sequence>&,
                const bool verbose = true);
}
#endif
