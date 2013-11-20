#include <iostream>
#include "gtest/gtest.h"
#include "gtr.hpp"
#include "sequence.hpp"

#include "Eigen/Core"

#include <Bpp/Phyl/Model.all>
#include <Bpp/Seq/Alphabet/DNA.h>

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

TEST(GTR, matches_bpp) {
    bpp::DNA dna;
    const Eigen::Vector4d pi(0.255068, 0.24877, 0.29809, 0.198071);
    const Eigen::Vector3d theta = gtr::piToTheta(pi);
    const std::vector<double> v{0.988033, 0.471959, 0.30081, 0.385086, 0.666584};
    bpp::GTR g(&dna, v[0], v[1], v[2], v[3], v[4], pi[0], pi[1], pi[2], pi[3]);
    gtr::GTRParameters p;
    p.theta = theta;
    for(size_t i = 0; i < v.size(); i++) {
        p.params[i] = v[i];
    }
    const Eigen::Matrix4d q = p.createQMatrix();

    for(size_t i = 0; i < 4; i++) {
        for(size_t j = 0; j < 4; j++) {
            EXPECT_NEAR(g.Qij(i, j), q(i, j), 1e-3);
        }
    }

    const gtr::GTRModel m = p.createModel();

    auto bpp_eval = g.getEigenValues();
    auto eig_eval = m.decomp.eigenvalues().real();
    for(size_t i = 0; i < 4; i++) {
        EXPECT_NEAR(bpp_eval[i], eig_eval[i], 1e-3);
    }
}
