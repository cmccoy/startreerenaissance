#include <iostream>
#include "gtest/gtest.h"
#include "star_optim.hpp"
#include "sequence.hpp"

#include "Eigen/Core"

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Sequence.h>

TEST(GTR, simple_jc) {
    bpp::DNA dna;
    bpp::GTR model(&dna);
    bpp::ConstantDistribution rates(1.0);

    std::vector<Sequence> v { Sequence() };
    Sequence& s = v.front();

    Eigen::Matrix4d& m = s.substitutions;
    m.fill(0.0);
    m.diagonal() = Eigen::Vector4d::Constant(1.0);

    s.distance = 0.02;

    const double ll = star_optim::starLikelihood(model, rates, v);
    const double expll = -5.62490959465585; // from bppml
    EXPECT_NEAR(expll, ll, 1e-5);
}

TEST(GTR, known_distance) {
    bpp::DNA dna;
    bpp::GTR model(&dna);
    bpp::ConstantDistribution rates(1.0);

    std::vector<Sequence> v { Sequence() };

    Sequence& s = v.front();
    s.substitutions <<
        94, 3, 2, 1,
        2, 95, 2, 1,
        2, 4, 89, 5,
        1, 3, 2, 94;
    s.distance = 0.02;

    const double expected = -702.603126357669; // From bppml
    const double actual = star_optim::starLikelihood(model, rates, v);
    ASSERT_NEAR(expected, actual, 1e-5);
}


TEST(GTR, gamma_variation) {
    bpp::DNA dna;
    bpp::GTR model(&dna);
    bpp::GammaDiscreteDistribution rates(4, 1.2, 1.0);
    model.setParameterValue("a", 0.5);
    model.setParameterValue("theta", 0.6);
    model.setParameterValue("theta1", 0.4);

    std::vector<Sequence> v { Sequence() };

    Sequence& s = v.front();
    s.substitutions <<
        94, 3, 2, 1,
        2, 95, 2, 1,
        2, 4, 89, 5,
        1, 3, 2, 94;
    s.distance = 0.02;

    // bppml input.sequence.file=test.fasta input.tree.file=test.tre output.tree.file=/dev/null model=GTR(a=0.5, theta=0.6) optimization=None rate_distribution=Gamma(n=4, alpha=1.2)
    const double expected = -714.208570235786; // From bppml
    const double actual = star_optim::starLikelihood(model, rates, v);
    ASSERT_NEAR(expected, actual, 1e-5);
}
