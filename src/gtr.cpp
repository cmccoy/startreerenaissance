#include "gtr.hpp"
#include "sequence.hpp"
#include <boost/math/tools/minima.hpp>
#include <nlopt.hpp>

using Eigen::Array4d;
using Eigen::Matrix4d;
using Eigen::Vector4d;


namespace gtr
{

/// Minimum improvement in LL over a round
const double IMPROVE_THRESH = 0.1;
const size_t MAX_ROUNDS = 200;
const size_t MAX_ITER = 300;
const double MIN_SUBS_PARAM = 1e-5,
             MAX_SUBS_PARAM = 20.0;
const size_t BIT_TOL = 50;

// GTRModel
Matrix4d GTRModel::createPMatrix(const double t) const
{
    const Matrix4d& v = decomp.eigenvectors().real();
    Vector4d lambda = (Array4d(decomp.eigenvalues().real()) * t).exp();
    return v * lambda.asDiagonal() * v.inverse();
}

double GTRModel::logLikelihood(const Sequence& s) const
{
    const Matrix4d p = createPMatrix(s.distance);
    //std::cout << "P(" << s.distance << ") = " << p << "\n\n";

    auto f = [](const double d) { return std::log(d); };
    const Matrix4d logP = p.unaryExpr(f);

    return logP.cwiseProduct(s.substitutions).sum();
}

// GTRParameters

GTRParameters::GTRParameters()
{
    params.fill(1);
    theta.fill(0.5);
};

/// \brief Create a rate matrix
///
/// This is parameterized as in Bio++: see
/// http://biopp.univ-montp2.fr/apidoc/bpp-phyl/html/classbpp_1_1GTR.html
Matrix4d GTRParameters::createQMatrix() const
{
    const Vector4d pi = createBaseFrequencies();
    const double a = params[0], b = params[1], c = params[2],
                 d = params[3], e = params[4], f = 1;
    Matrix4d S;
    S << 0, d, f, b,
         d, 0, e, a,
         f, e, 0, c,
         b, a, c, 0;

    const Vector4d diag((-b * pi[3] - pi[2] - d * pi[1]) / pi[0],
                        (-a * pi[3] - e * pi[2] - d * pi[0]) / pi[1],
                        (-c * pi[3] - e * pi[1] - pi[0]) / pi[2],
                        (-c * pi[2] - a * pi[1] - b * pi[0]) / pi[3]);
    S.diagonal() = diag;

    // Normalization
    const double p = 2 * (a * pi[1] * pi[3] + b * pi[0] * pi[3] + c * pi[2] * pi[3] + d * pi[0] * pi[1] + e * pi[1] * pi[2] + f * pi[0] * pi[2]);
    S /= p;


    //assert(S.rowwise().sum().isZero());
    const Matrix4d Q = S * pi.asDiagonal();
    // Sanity check: this is just scaling, can be approximate.
    //if(std::abs(Q.diagonal().cwiseProduct(pi).sum() + 1) > 0.1) {
        //std::cerr << "S=" << S << '\n';
        //std::cerr << "Q=" << Q << '\n';
        //std::cerr << "theta=" << theta.transpose() << '\n';
        //std::cerr << "pi=" << pi.transpose() << '\n';
        //std::cerr << "sum(diag(Q)*pi) = " << Q.diagonal().cwiseProduct(pi).sum() << '\n';
        //assert(0 && "Constraint violated.");
    //}

    return Q;
};

Eigen::Vector4d thetaToPi(const Eigen::Vector3d& theta)
{
    Vector4d r;
    r << theta[1] * (1 - theta[0]), (1 - theta[2]) * theta[0], theta[2] * theta[0], (1 - theta[1]) * (1 - theta[0]);
    return std::move(r);
}

Eigen::Vector3d piToTheta(Eigen::Vector4d pi)
{
    pi /= pi.sum();
    Eigen::Vector3d result;
    result << pi[1] + pi[2], pi[0] / (pi[0] + pi[3]), pi[2] / (pi[1] + pi[2]);
    return result;
}

Vector4d GTRParameters::createBaseFrequencies() const
{
    return thetaToPi(theta);
}

GTRModel GTRParameters::createModel() const
{
    return GTRModel(createQMatrix());
}

double& GTRParameters::parameter(size_t index) {
    assert(index < numberOfParameters() && "Invalid index");
    const size_t offset = 5;
    if(index < offset) {
        return params[index];
    } else {
        return theta[index - offset];
    }
}

double GTRParameters::parameter(size_t index) const {
    assert(index < numberOfParameters() && "Invalid index");
    const size_t offset = 5;
    if(index < offset) {
        return params[index];
    } else {
        return theta[index - offset];
    }
}

std::vector<double> GTRParameters::toVector() const
{
    std::vector<double> r(numberOfParameters());
    for(size_t i = 0; i < numberOfParameters(); i++)
        r[i] = parameter(i);
    return r;
}

void GTRParameters::ofVector(const std::vector<double>& v)
{
    assert(v.size() == numberOfParameters() && "Parameter size mismatch");
    for(size_t i = 0; i < numberOfParameters(); i++)
        parameter(i) = v[i];
}

// Functions
double starLikelihood(const GTRModel& model,
                      const std::vector<Sequence>& sequences)
{
    double result = 0.0;

    #pragma omp parallel for reduction(+:result)
    for(size_t i = 0; i < sequences.size(); i++) {
        result += model.logLikelihood(sequences[i]);
    }
    return result;
}

void estimateBranchLengths(const GTRModel& model,
                           std::vector<Sequence>& sequences)
{
    #pragma omp parallel for
    for(size_t i = 0; i < sequences.size(); i++) {
        Sequence& s = sequences[i];
        auto f = [&model, &s](const double d) {
            s.distance = d;
            const double result = -model.logLikelihood(s);
            return result;
        };
        boost::uintmax_t max_iter = 100;
        std::pair<double, double> res =
            boost::math::tools::brent_find_minima(f, 1e-9, 1.0, 50, max_iter);
        s.distance = res.first;
    }
}

Eigen::Vector3d countBaseFrequencies(const std::vector<Sequence>& sequences)
{
    Eigen::Vector4d counts;
    counts.fill(1); // Pseudocount

    for(const Sequence & s : sequences)
        counts += s.substitutions.colwise().sum();

    Eigen::Vector4d pi = counts / counts.sum();

    Eigen::Vector3d result;
    result << pi[1] + pi[2], pi[0] / (pi[0] + pi[3]), pi[2] / (pi[1] + pi[2]);
    assert((thetaToPi(result) - pi).isZero() && "pi to theta conversion failed");
    return result;
}

void empiricalModel(const std::vector<Sequence>& sequences,
                    gtr::GTRParameters& model)
{
    model.theta = countBaseFrequencies(sequences);

    Matrix4d result;
    result.fill(0);

    for(const Sequence & s : sequences) {
        result += s.substitutions;
    }

    model.params << result(1, 3), result(0, 3), result(2, 3),
                    result(0, 1), result(1, 2);
    model.params /= model.params.sum();
}

double optimizeParameter(const std::vector<Sequence>& sequences,
                         const size_t index,
                         gtr::GTRParameters& params,
                         const double min_value=MIN_SUBS_PARAM,
                         const double max_value=MAX_SUBS_PARAM)
{
    double& parameter = params.parameter(index);
    auto f = [&sequences, &parameter, &params](const double d) {
        parameter = d;
        const GTRModel model = params.createModel();
        return -starLikelihood(model, sequences);
    };

    boost::uintmax_t max_iter = MAX_ITER;
    std::pair<double, double> result =
        boost::math::tools::brent_find_minima(f, min_value, max_value, BIT_TOL, max_iter);

    parameter = result.first;

    return -result.second;
}

struct NlOptParams
{
    std::vector<Sequence>* sequences;
    GTRParameters* params;
};

double nlLogLike(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
    assert(grad.empty());

    NlOptParams* params = reinterpret_cast<NlOptParams*>(data);

    params->params->ofVector(x);

    return starLikelihood(params->params->createModel(), *params->sequences);
}

// TODO: tolerance, no magic numbers, no printing.
void optimize(gtr::GTRParameters& params,
              std::vector<Sequence>& sequences,
              bool verbose)
{
    double lastLogLike = starLikelihood(params.createModel(), sequences);
    const size_t nParam = params.numberOfParameters();

    for(size_t iter = 0; iter < MAX_ROUNDS; iter++) {
        bool anyImproved = false;

        // First, substitution model
        nlopt::opt opt(nlopt::LN_BOBYQA, nParam);
        NlOptParams optParams { &sequences, &params };
        opt.set_max_objective(nlLogLike, &optParams);

        std::vector<double> lowerBounds(nParam, MIN_SUBS_PARAM);
        std::vector<double> upperBounds(nParam, MAX_SUBS_PARAM);
        for(size_t i = 5; i < nParam; i++) {
            lowerBounds[i] = 0.01;
            upperBounds[i] = 0.99;
        }
        opt.set_lower_bounds(lowerBounds);
        opt.set_upper_bounds(upperBounds);
        opt.set_ftol_abs(IMPROVE_THRESH);
        opt.set_maxeval(MAX_ITER);

        std::vector<double> x = params.toVector();

        double logLike;
        const int nlOptResult = opt.optimize(x, logLike);
        //std::clog << "Optimization finished with " << nlOptResult << '\n';
        //if(nlOptResult != nlopt::SUCCESS)

        params.ofVector(x);

        if(verbose) {
            std::clog << "iteration " << iter << ": " << lastLogLike << " ->\t" << logLike << '\t' << logLike - lastLogLike << '\n';
            std::clog << "p=" << params.params.transpose() << '\t' << "theta=" << params.theta.transpose() << '\n';
            std::clog.flush();
        }

        if(std::abs(logLike - lastLogLike) > IMPROVE_THRESH)
            anyImproved = true;
        lastLogLike = logLike;


        // then, branch lengths
        estimateBranchLengths(params.createModel(),
                              sequences);
        logLike = starLikelihood(params.createModel(), sequences);

        if(verbose) {
            std::clog << "iteration " << iter << " (branch lengths): " << lastLogLike << " ->\t" << logLike << '\t' << logLike - lastLogLike << '\n';
            std::clog.flush();
        }

        if(std::abs(logLike - lastLogLike) > IMPROVE_THRESH)
            anyImproved = true;
        lastLogLike = logLike;


        if(!anyImproved)
            break;
    }
}

}
