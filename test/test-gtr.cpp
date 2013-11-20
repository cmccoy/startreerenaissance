#include <iostream>
#include "gtest/gtest.h"
#include "gtr.hpp"
#include "sequence.hpp"

#include "Eigen/Core"

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Sequence.h>

TEST(GTR, simple_jc) {
    gtr::GTRParameters params;
    Sequence s;
    Eigen::Matrix4d& m = s.substitutions;
    m.fill(0.0);
    m.diagonal() = Eigen::Vector4d::Constant(1.0);

    s.distance = 0.02;

    const double ll = params.createModel().logLikelihood(s);
    const double expll = -0.07973215; // From by hand calculation in R
    EXPECT_NEAR(expll, ll, 1e-5);
}

TEST(GTR, known_distance) {
    Sequence s;
    s.substitutions << 
        94, 3, 2, 1,
        2, 95, 2, 1,
        2, 4, 89, 5,
        1, 3, 2, 94;
    s.distance = 0.02;

    gtr::GTRParameters parameters;
    const double initLl = parameters.createModel().logLikelihood(s);
    ASSERT_NEAR(-148.0854, initLl, 0.5); // From by-hand calculation in R - see test-calcs.R
}

TEST(GTR, roundtrip) {
    Eigen::Vector4d pi;
    pi << 0.244006, 0.224878, 0.297143, 0.233973;

    const Eigen::Vector3d theta = gtr::piToTheta(pi);

    const Eigen::Vector4d pi_conv = gtr::thetaToPi(theta);

    for(size_t i = 0; i < 4; i++) {
        EXPECT_NEAR(pi[i], pi_conv[i], 1e-4);
    }
}

class BppCompare : public ::testing::Test 
{
protected:
    BppCompare() :
        g(&dna) {};
    virtual void SetUp() {
        const Eigen::Vector4d pi(0.255068, 0.24877, 0.29809, 0.198071);
        const Eigen::Vector3d theta = gtr::piToTheta(pi);
        const std::vector<double> v{0.988033, 0.471959, 0.30081, 0.385086, 0.666584};
        g = bpp::GTR(&dna, v[0], v[1], v[2], v[3], v[4], pi[0], pi[1], pi[2], pi[3]);
        p.theta = theta;
        for(size_t i = 0; i < v.size(); i++) {
            p.params[i] = v[i];
        }
    }

    // virtual void TearDown() {}

    bpp::DNA dna;
    bpp::GTR g;
    gtr::GTRParameters p;
};

void check_against_bpp(double a, double b, double c, double d, double e, double p1, double p2, double p3, double p4, double time)
{

    const double tol = 1e-4;
    bpp::DNA dna;
    const Eigen::Vector4d pi(p1, p2, p3, p4);
    const Eigen::Vector3d theta = gtr::piToTheta(pi);
    bpp::GTR g(&dna, a, b, c, d, e, p1, p2, p3, p4);
    gtr::GTRParameters p;

    p.params << a, b, c, d, e;
    p.theta = theta;

    std::vector<std::string> names {"a", "b", "c", "d", "e", "theta", "theta1", "theta2"};

    for(size_t i = 0; i < names.size(); i++) {
        EXPECT_NEAR(g.getParameterValue(names[i]), p.parameter(i), tol) << "Parameter " << names[i] << " not equal";
    }

    const Eigen::Matrix4d q = p.createQMatrix();

    const size_t dim = 4;
    for(size_t i = 0; i < dim; i++) {
        for(size_t j = 0; j < dim; j++) {
            ASSERT_NEAR(g.Qij(i, j), q(i, j), tol) << "Q(" << i << ", " << j << ")";
        }
    }

    const auto bpp_eval = g.getEigenValues();
    const gtr::GTRModel m = p.createModel();
    const auto eig_eval = m.decomp.eigenvalues().real();

    for(int i = 0; i < eig_eval.size(); i++) {
        ASSERT_NEAR(bpp_eval[i], eig_eval[i], tol) << "Eigenvalue " << i;
    }

    const Eigen::Matrix4d p_eig = m.createPMatrix(time);
    const bpp::Matrix<double>& p_bpp = g.getPij_t(time);
    ASSERT_EQ(dim, p_bpp.getNumberOfRows());
    ASSERT_EQ(dim, p_bpp.getNumberOfColumns());
    for(size_t i = 0; i < dim; i++) {
        const std::vector<double> r = p_bpp.row(i);
        for(size_t j = 0; j < dim; j++) {
            ASSERT_NEAR(r[j], p_eig(i, j), tol) << "P(" << i << "," << j << ")";
        }
    }
}

TEST(BppCompareDecomp, decomposition_matches_bpp1) {
    const Eigen::Vector4d pi(0.255068, 0.24877, 0.29809, 0.198071);
    const std::vector<double> v{0.988033, 0.471959, 0.30081, 0.385086, 0.666584};

    for(double t : {0.1, 0.01, 0.4})
        check_against_bpp(0.988033, 0.471959, 0.30081, 0.385086, 0.666584, 0.255068, 0.24877, 0.29809, 0.198071, t);
}

TEST_F(BppCompare, distance_estimation) {
    Sequence s;
    s.substitutions << 
        94, 3, 2, 1,
        2, 95, 2, 1,
        2, 4, 89, 5,
        1, 3, 2, 94;
    s.distance = 0.02;

    std::vector<std::string> sequences(2);
    const std::vector<char> bases {'A', 'C', 'G', 'T'};
    for(size_t i = 0; i < 4; i++) {
        for(size_t j = 0; j < 4; j++) {
            const size_t n = s.substitutions(i, j);
            for(size_t k = 0; k < n; k++) {
                sequences[0] += bases[i];
                sequences[1] += bases[j];
            }
        }
    }

    std::vector<Sequence> seqs { s };
    const gtr::GTRModel model = p.createModel();
    gtr::estimateBranchLengths(model, seqs);

    bpp::VectorSiteContainer sites(&dna);
    sites.addSequence(bpp::BasicSequence("ref", sequences[0], &dna));
    sites.addSequence(bpp::BasicSequence("query", sequences[1], &dna));

    bpp::ConstantDistribution rate(1.0);
    bpp::DistanceEstimation est(&g, &rate, &sites, 0);
    est.computeMatrix();

    EXPECT_NEAR(seqs[0].distance, est.getMatrix()->operator()("ref", "query"), 1e-3);

    // Estimate all parameters
    //gtr::optimize(p, seqs, false);
    //std::cerr << "log(L) = " << gtr::starLikelihood(p.createModel(), seqs) << '\n';
    //est.setAdditionalParameters(g.getParameters());
    //est.computeMatrix();

    //const std::vector<std::string> names { "a", "b", "c", "d", "e", "theta", "theta1", "theta2" };
    //const std::vector<double> gtr { p.params[0], p.params[1], p.params[2], p.params[3],
                                    //p.params[4], p.theta[0], p.theta[1], p.theta[2] };
    //for(size_t i = 0; i < names.size(); i++) {
        //EXPECT_NEAR(g.getParameterValue(names[i]), gtr[i], 0.05) << names[i];
    //}
    //EXPECT_NEAR(seqs[0].distance, est.getMatrix()->operator()("ref", "query"), 1e-3);
}
