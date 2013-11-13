#include <iostream>
#include "gtest/gtest.h"
#include "gtr.hpp"
#include "sequence.hpp"

#include "Eigen/Core"

TEST(simple_jc, simple_jc) {
    gtr::GTRParameters params;
    params.params.fill(1.0);
    Sequence s;
    Eigen::Matrix4d& m = s.substitutions;
    m.fill(0.0);
    m.diagonal() = Eigen::Vector4d::Constant(1.0);

    s.distance = 0.01;

    const double ll = params.buildModel().logLikelihood(s);
    const double expll = std::log(0.9705921) * 4;
    EXPECT_NEAR(expll, ll, 1e-5);
}
