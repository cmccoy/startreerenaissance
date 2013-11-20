#ifndef GTR_GTR_H
#define GTR_GTR_H

#include <vector>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

struct Sequence;

namespace gtr
{

using Vector5d = Eigen::Matrix<double, 5, 1>;
using EigenDecomp = Eigen::EigenSolver<Eigen::Matrix4d>;

struct GTRModel {
    Eigen::Matrix4d model;
    EigenDecomp decomp;

    GTRModel(Eigen::Matrix4d model) : model(model), decomp(model, true)
    { };

    /// \brief Get the P matrix for time t
    Eigen::Matrix4d createPMatrix(const double t) const;
    double logLikelihood(const Sequence& s) const;
};

struct GTRParameters {
    GTRParameters();

    Eigen::Matrix4d createQMatrix() const;
    Eigen::Vector4d createBaseFrequencies() const;
    GTRModel createModel() const;

    // 5 parameters of the GTR
    Vector5d params;
    // Length-3 base-frequency vector
    // as http://biopp.univ-montp2.fr/apidoc/bpp-phyl/html/classbpp_1_1GTR.html
    Eigen::Vector3d theta;

    inline size_t numberOfParameters() const { return 8; }
    double& parameter(size_t index);
    double parameter(size_t index) const;

    std::vector<double> toVector() const;
    void ofVector(const std::vector<double>&);

};

/// \brief calculate the star-tree likelihood
double starLikelihood(const GTRModel&,
                      const std::vector<Sequence>&);

/// \brief estimate branch lengths
void estimateBranchLengths(const GTRModel&,
                           std::vector<Sequence>&);

/// \brief Guess a model by counting
void empiricalModel(const std::vector<Sequence>&,
                    gtr::GTRParameters&);

/// \brief Optimize the GTR model & branch lengths for a collection of sequences
void optimize(gtr::GTRParameters&, std::vector<Sequence>&, bool verbose=true);

/// \brief Convert from three-parameter to 4-parameter
Eigen::Vector4d thetaToPi(const Eigen::Vector3d&);
Eigen::Vector3d piToTheta(Eigen::Vector4d pi);

}
#endif
