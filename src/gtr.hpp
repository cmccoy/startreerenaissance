#ifndef GTR_GTR_H
#define GTR_GTR_H

#include <string>
#include <vector>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

struct Sequence;

namespace gtr {

typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::EigenSolver<Eigen::Matrix4d> EigenDecomp;

struct GTRModel
{
    Eigen::Matrix4d model;
    EigenDecomp decomp;

    GTRModel(Eigen::Matrix4d model) : model(model), decomp(model, true)
    { };

    /// \brief Get the P matrix for time t
    Eigen::Matrix4d buildPMatrix(const double t) const;
    double logLikelihood(const Sequence& s) const;
};

struct GTRParameters
{
    GTRParameters();

    Eigen::Matrix4d buildQMatrix() const;
    GTRModel buildModel() const;

    // 6 parameters of the GTR
    Vector6d params;
    // base frequencies
    Eigen::Vector4d pi;
};

}

#endif
