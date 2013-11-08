#ifndef GTR_GTR_H
#define GTR_GTR_H

#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace gtr {

using Eigen::Matrix4d;
using Eigen::Vector4d;
using Eigen::Array4d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::EigenSolver<Matrix4d> EigenDecomp;

struct Sequence
{
    Matrix4d transitions;
    double distance;
};

struct GTRModel
{
    Matrix4d model;
    EigenDecomp decomp;

    GTRModel(Matrix4d model) : model(model), decomp(model, true)
    { };

    /// \brief Get the P matrix for time t
    Matrix4d buildPMatrix(const double t) const;
    double logLikelihood(const Sequence& s) const;
};

struct GTRParameters
{
    GTRParameters();

    Matrix4d buildQMatrix() const;
    GTRModel buildModel() const;

    // 6 parameters of the GTR
    Vector6d params;
    // base frequencies
    Vector4d pi;
};

}

#endif
