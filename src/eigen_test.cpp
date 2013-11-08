#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/SVD>


using Eigen::Matrix4d;
using Eigen::Vector4d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::JacobiSVD<Matrix4d> SubsMatrixDecomp;

struct Sequence
{
    Matrix4d transitions;
    double distance;
};

struct GTRModel
{
    GTRModel(Matrix4d model) : model(model), decomp(model)
    {};
    Matrix4d model;
    SubsMatrixDecomp decomp;
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

int main()
{
    GTRParameters m;
    GTRModel model = m.buildModel();

    std::cout << "singular values:\n";
    std::cout << model.decomp.singularValues() << '\n';

    std::cout << m.params << '\n';
    std::cout << m.pi << '\n';
    std::cout << m.buildQMatrix() << '\n';
    return 0;
}
