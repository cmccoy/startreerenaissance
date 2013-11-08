#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/SVD>


using Eigen::Matrix4d;
using Eigen::Vector4d;
using Eigen::Array4d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::JacobiSVD<Matrix4d> SubsMatrixDecomp;


struct GTRModel
{
    GTRModel(Matrix4d model) : model(model), decomp(model)
    {};
    Matrix4d model;
    SubsMatrixDecomp decomp;

    Matrix4d buildPMatrix(const double t) const
    {
        const Matrix4d& u = decomp.matrixU();
        const Matrix4d& v = decomp.matrixV();
        Vector4d lambda = (Array4d(decomp.singularValues()) * t).exp();
        return u * Eigen::DiagonalMatrix<double, 4, 4>(lambda) * v;
    }
};

struct GTRParameters
{
    GTRParameters()
    {
        params.fill(1);
        pi.fill(0.25);
    };

    Matrix4d buildQMatrix() const
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

    GTRModel buildModel() const
    {
        return GTRModel(buildQMatrix());
    }

    // 6 parameters of the GTR
    Vector6d params;
    // base frequencies
    Vector4d pi;
};

struct Sequence
{
    Matrix4d transitions;
    double distance;

    double logLikelihood(const GTRModel& model) const
    {
        const Matrix4d u = model.decomp.matrixU();
        const Matrix4d v = model.decomp.matrixV();
        Array4d lambda = model.decomp.singularValues();
        lambda = lambda * distance;
        lambda = lambda.exp();
        Vector4d lambdaExp = lambda;
        const Matrix4d p = u * Eigen::DiagonalMatrix<double, 4, 4>(lambdaExp) * v * transitions;
        double result = 0;
        for(size_t i = 0; i < 4; i++) {
            for(size_t j = 0; j < 4; j++) {
                const double d = p(i, j);
                if(d > 0)
                    result += std::log(d);
            }
        }

        return result;
    };
};


int main()
{
    GTRParameters m;
    GTRModel model = m.buildModel();

    std::cout << "singular values:\n";
    std::cout << model.decomp.singularValues() << '\n';

    std::cout << m.params << '\n';
    std::cout << m.pi << '\n';
    std::cout << m.buildQMatrix() << '\n';
    std::cout << "Distance = 1\n" << model.buildPMatrix(1.0) << '\n';
    return 0;
}
