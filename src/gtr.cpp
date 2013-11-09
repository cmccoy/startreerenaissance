#include "gtr.hpp"

namespace gtr {

// GTRModel
Matrix4d GTRModel::buildPMatrix(const double t) const
{
    const Matrix4d& v = decomp.eigenvectors().real();
    Vector4d lambda = (Array4d(decomp.eigenvalues().real()) * t).exp();
    return v * Eigen::DiagonalMatrix<double, 4, 4>(lambda) * v.inverse();
}

double GTRModel::logLikelihood(const Sequence& s) const
{
    const Matrix4d p = buildPMatrix(s.distance);

    auto f = [](const double d) { return std::log(d); };
    const Matrix4d log_p = p.unaryExpr(f);

    return log_p.cwiseProduct(s.transitions).sum();
}

// GTRParameters

GTRParameters::GTRParameters()
{
    params.fill(1);
    pi.fill(0.25);
};

Matrix4d GTRParameters::buildQMatrix() const
{
    const double pi1 = pi[0], pi2 = pi[1], pi3 = pi[2], pi4 = pi[3];
    const double x1 = params[0], x2 = params[1], x3 = params[2],
          x4 = params[3], x5 = params[4], x6 = params[5];
    Matrix4d result;
    result << - (x1 + x2 + x3), (pi1 * x1 / pi2), (pi1 * x2 / pi3), pi1 * x3 / pi4,
           x1, -(pi1 * x1 / pi2 + x4 + x5), (pi2 * x4 / pi3), (pi2 * x5 / pi4),
           x2, x4, -(pi1 * x2 / pi3 + pi2 * x4 / pi3 + x6), pi3 * x6 / pi4,
           x3, x5, x6, -(pi1 * x3 / pi4 + pi2 * x5 / pi4 + pi3 * x6 / pi4);
    return result;
};

GTRModel GTRParameters::buildModel() const
{
    return GTRModel(buildQMatrix());
}

}
