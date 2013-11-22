#include <iostream>
#include "gtest/gtest.h"
#include "star_optim.hpp"
#include "sequence.hpp"

#include "Eigen/Core"

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Sequence.h>

bpp::DNA DNA;

bpp::VectorSiteContainer createSites(const Sequence& sequence)
{
    std::vector<std::string> names {"ref", "qry" };
    std::vector<std::string> seqs(2);
    const std::string bases = "ACGT";

    for(size_t i = 0; i < bases.size(); i++) {
        for(size_t j = 0; j < bases.size(); j++) {
            int count = static_cast<int>(sequence.substitutions(i, j));
            for(int k = 0; k < count; k++) {
                seqs[0] += bases[i];
                seqs[1] += bases[j];
            }
        }
    }

    bpp::VectorSiteContainer result(&DNA);
    for(int i = 0; i < 2; i++)
        result.addSequence(bpp::BasicSequence(names[i], seqs[i], &DNA));

    return result;
}

void checkAgainstBpp(const std::vector<Sequence>& sequences,
                       bpp::SubstitutionModel& model,
                       bpp::DiscreteDistribution& rates)
{
    using namespace bpp;
    ASSERT_EQ(1, sequences.size());
    const double starLL = star_optim::starLikelihood(model, rates, sequences);

    VectorSiteContainer sites = createSites(sequences[0]);

    Node *root = new Node(0),
         *c1 = new Node(1, sites.getSequence(0).getName()),
         *c2 = new Node(2, sites.getSequence(1).getName());
    root->addSon(c1);
    root->addSon(c2);
    c1->setDistanceToFather(sequences[0].distance / 2);
    c2->setDistanceToFather(sequences[0].distance / 2);
    TreeTemplate<Node> tree(root);

    DRHomogeneousTreeLikelihood calc(tree, sites, &model, &rates, true, false);
    calc.initialize();
    calc.computeTreeLikelihood();

    const double bppLL = calc.getLogLikelihood();
    EXPECT_NEAR(bppLL, starLL, 1e-3) << "Likelihood calculations do not match.";
}

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

    checkAgainstBpp(v, model, rates);

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

    checkAgainstBpp(v, model, rates);
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

    checkAgainstBpp(v, model, rates);
}
