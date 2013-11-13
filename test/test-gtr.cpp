#include <iostream>
#include "gtest/gtest.h"
#include "gtr.hpp"
#include "sequence.hpp"

#include "Eigen/Core"

TEST(simple_jc, simple_jc) {
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

TEST(known_distance, known_distance) {
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
