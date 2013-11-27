#include "star_optim.hpp"
#include "sequence.hpp"

#ifdef HAVE_OMP
#include <omp.h>
#endif

#include "libhmsbeagle/beagle.h"

#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/tools/minima.hpp>

#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>

#include <string>
#include <stdexcept>

#include <nlopt.hpp>

namespace star_optim
{

/// Minimum improvement in LL over a round
const double IMPROVE_THRESH = 0.1;
const size_t MAX_ROUNDS = 30;
const size_t MAX_ITER = 300;
const double MIN_SUBS_PARAM = 1e-5,
             MAX_SUBS_PARAM = 20.0;
const size_t BIT_TOL = 50;

void beagleCheck(const int value, const std::string& details = "")
{
    std::string s;
    switch(value) {
        case BEAGLE_ERROR_GENERAL: s = "BEAGLE_ERROR_GENERAL"; break;
        case BEAGLE_ERROR_OUT_OF_MEMORY: s = "BEAGLE_ERROR_OUT_OF_MEMORY"; break;
        case BEAGLE_ERROR_UNIDENTIFIED_EXCEPTION: s = "BEAGLE_ERROR_UNIDENTIFIED_EXCEPTION"; break;
        case BEAGLE_ERROR_UNINITIALIZED_INSTANCE: s = "BEAGLE_ERROR_UNINITIALIZED_INSTANCE"; break;
        case BEAGLE_ERROR_OUT_OF_RANGE: s = "BEAGLE_ERROR_OUT_OF_RANGE"; break;
        case BEAGLE_ERROR_NO_RESOURCE: s = "BEAGLE_ERROR_NO_RESOURCE"; break;
        case BEAGLE_ERROR_NO_IMPLEMENTATION: s = "BEAGLE_ERROR_NO_IMPLEMENTATION"; break;
        case BEAGLE_ERROR_FLOATING_POINT: s = "BEAGLE_ERROR_FLOATING_POINT"; break;
        default: break;
    }
    if(!s.empty())
        throw std::runtime_error(s + " " + details);
}



/// Copy the contents of vec into arr
/// \param arr destination array, with length at least vec.size()
/// \param vec Vector to copy from
void blitVectorToArray(double* arr, const std::vector<double>& vec)
{
    for(std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it)
        *arr++ = *it;
}

/// Copy the contents of matrix into arr, in row-major order
/// \param arr destination array, with length at least nrows x ncols in length
/// \param matrix Vector to copy from
void blitMatrixToArray(double* arr, const bpp::Matrix<double>& matrix)
{
    const int cols = matrix.getNumberOfColumns(), rows = matrix.getNumberOfRows();
    for(int i = 0; i < rows; ++i) {
        blitVectorToArray(arr, matrix.row(i));
        arr += cols;
    }
}

double pairLogLikelihood(const int beagleInstance, const Sequence& sequence, const int nStates)
{
    std::vector<double> patternWeights;
    patternWeights.reserve(nStates * nStates);
    for(int i = 0; i < nStates; i++) {
        for(int j = 0; j < nStates; j++) {
            patternWeights.push_back(sequence.substitutions(i, j));
        }
    }

    beagleCheck(beagleSetPatternWeights(beagleInstance, patternWeights.data()));

    std::vector<int> nodeIndices { 0, 1 };
    std::vector<double> edgeLengths { sequence.distance, 0 };
    const int rootIndex = 2;

    beagleCheck(beagleUpdateTransitionMatrices(beagleInstance,
                0,
                nodeIndices.data(),
                nullptr,
                nullptr,
                edgeLengths.data(),
                nodeIndices.size()));

    BeagleOperation op {rootIndex, BEAGLE_OP_NONE, BEAGLE_OP_NONE,
                        nodeIndices[0], nodeIndices[0],
                        nodeIndices[1], nodeIndices[1]
                       };

    beagleCheck(beagleUpdatePartials(beagleInstance,      // instance
                                     &op,
                                     1,  // # of ops
                                     rootIndex), // cumulative scaling index
                "updatePartials");

    int scaleIndex = op.destinationPartials;
    beagleAccumulateScaleFactors(beagleInstance, &scaleIndex, 1, BEAGLE_OP_NONE);
    double logLike = 1;
    const int weightIndex = 0, freqIndex = 0;
    int returnCode = beagleCalculateRootLogLikelihoods(beagleInstance,               // instance
                     &rootIndex,
                     &weightIndex,
                     &freqIndex,
                     &rootIndex,
                     1,
                     &logLike);
    beagleCheck(returnCode, "rootLogLike");
    return logLike;
}

int createBeagleInstance(const bpp::SubstitutionModel& model,
                         const bpp::DiscreteDistribution& rates)
{
    const int nStates = model.getNumberOfStates();
    const int nRates = rates.getNumberOfCategories();
    const int nTips = 2;
    const int nBuffers = 3;
    BeagleInstanceDetails deets;
    const int instance = beagleCreateInstance(nTips,                            /**< Number of tip data elements (input) */
                         nBuffers,       /**< Number of partials buffers to create (input) */
                         nTips,                    /**< Number of compact state representation buffers to create (input) */
                         nStates,           /**< Number of states in the continuous-time Markov chain (input) */
                         nStates * nStates,            /**< Number of site patterns to be handled by the instance (input) */
                         1,                    /**< Number of rate matrix eigen-decomposition buffers to allocate (input) */
                         nBuffers,                    /**< Number of rate matrix buffers (input) */
                         nRates,            /**< Number of rate categories (input) */
                         nBuffers,                       /**< Number of scaling buffers */
                         nullptr,                     /**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
                         0,                        /**< Length of resourceList list (input) */
                         BEAGLE_FLAG_VECTOR_SSE | BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_SCALING_AUTO, // Bit-flags indicating
                         0,                /**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
                         &deets);

    beagleCheck(instance);
    // Fill rates
    beagleSetCategoryRates(instance, rates.getCategories().data());
    beagleSetCategoryWeights(instance, 0, rates.getProbabilities().data());

    // And states
    std::vector<int> ref(nStates * nStates), qry(nStates * nStates);
    for(int i = 0; i < nStates; i++) {
        for(int j = 0; j < nStates; j++) {
            ref[nStates * i + j] = i;
            qry[nStates * i + j] = j;
        }
    }
    beagleSetTipStates(instance, 0, ref.data());
    beagleSetTipStates(instance, 1, qry.data());

    // And eigen decomposition
    std::vector<double> evec(nStates * nStates),
        ivec(nStates * nStates),
        eval(nStates);
    blitMatrixToArray(evec.data(), model.getColumnRightEigenVectors());
    blitMatrixToArray(ivec.data(), model.getRowLeftEigenVectors());
    blitVectorToArray(eval.data(), model.getEigenValues());
    beagleSetEigenDecomposition(instance, 0, evec.data(), ivec.data(), eval.data());

    // And state frequencies
    beagleSetStateFrequencies(instance, 0, model.getFrequencies().data());

    return instance;
}

double starLikelihood(const bpp::SubstitutionModel& model,
                      const bpp::DiscreteDistribution& rates,
                      const std::vector<Sequence>& sequences)
{
    const bpp::ParameterList p = model.getIndependentParameters();
    //for(std::size_t i = 0; i < p.size(); i++) {
    //std::clog << [>p[i].getName() << "=" <<<] p[i].getValue() << "\t";
    //}

    const int nStates = model.getNumberOfStates();
    int instance = -1;

    std::vector<int> beagleInstanceIds;
#ifdef HAVE_OMP
    omp_lock_t writelock;
    omp_init_lock(&writelock);
#endif

    double result = 0.0;
    #pragma omp parallel for reduction(+:result) firstprivate(instance)
    for(size_t i = 0; i < sequences.size(); i++) {
        if(instance <= 0) {
#ifdef HAVE_OMP
            omp_set_lock(&writelock);
#endif
            instance = createBeagleInstance(model, rates);
            beagleInstanceIds.push_back(instance);
#ifdef HAVE_OMP
            omp_unset_lock(&writelock);
#endif
        }
        result += pairLogLikelihood(instance, sequences[i], nStates);
    }
#ifdef HAVE_OMP
    omp_destroy_lock(&writelock);
#endif

    for(const int i : beagleInstanceIds)
        beagleFinalizeInstance(i);

    double prior = 0.0;
    assert(model.hasParameter("kappa"));
    if(model.hasParameter("kappa")) {
        boost::math::lognormal_distribution<double> distn(1, 1.25);
        prior += std::log(boost::math::pdf(distn, model.getParameterValue("kappa")));
    }

    return result + prior;
}

/// \brief estimate branch lengths
void estimateBranchLengths(const bpp::SubstitutionModel& model,
                           const bpp::DiscreteDistribution& rates,
                           std::vector<Sequence>& sequences)
{
    int instance = -1;

    std::vector<int> beagleInstanceIds;

    const size_t nStates = model.getNumberOfStates();

#ifdef HAVE_OMP
    omp_lock_t writelock;
    omp_init_lock(&writelock);
    #pragma omp parallel for firstprivate(instance)
#endif
    for(size_t i = 0; i < sequences.size(); i++) {
        if(instance <= 0) {
#ifdef HAVE_OMP
            omp_set_lock(&writelock);
#endif
            instance = createBeagleInstance(model, rates);
            beagleInstanceIds.push_back(instance);
#ifdef HAVE_OMP
            omp_unset_lock(&writelock);
#endif
        }

        Sequence& s = sequences[i];

        auto f = [instance, &s, nStates](const double d) {
            assert(!std::isnan(d) && "NaN distance?");
            s.distance = d;
            const double result = pairLogLikelihood(instance, s, nStates);
            return -result;
        };

        boost::uintmax_t max_iter = 100;
        std::pair<double, double> res =
            boost::math::tools::brent_find_minima(f, 1e-6, 0.8, 50, max_iter);
        assert(!std::isnan(res.first) && "NaN distance?");
        s.distance = res.first;
    }
#ifdef HAVE_OMP
    omp_destroy_lock(&writelock);
#endif

    for(const int i : beagleInstanceIds)
        beagleFinalizeInstance(i);
}

struct NlOptParams {
    std::vector<Sequence>* sequences;
    bpp::ParameterList* paramList;
    bpp::SubstitutionModel* model;
    bpp::DiscreteDistribution* rates;
};

double nlLogLike(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    assert(grad.empty() && "Expected no derivative");

    NlOptParams* params = reinterpret_cast<NlOptParams*>(data);

    for(size_t i = 0; i < x.size(); i++) {
        bpp::Parameter& p = params->paramList->operator[](i);
        p.setValue(x[i]);
    }

    params->model->matchParametersValues(*params->paramList);
    params->rates->matchParametersValues(*params->paramList);

    return starLikelihood(*params->model, *params->rates, *params->sequences);
}

/// \brief Optimize the model & branch lengths distribution for a collection of sequences
size_t optimize(bpp::SubstitutionModel& model,
                bpp::DiscreteDistribution& rates,
                std::vector<Sequence>& sequences,
                const bool verbose)
{
    bpp::ParameterList params = model.getIndependentParameters();
    // TODO: this is a crude hack to handle gamma distributed rates, only
    if(rates.hasParameter("alpha"))
        params.addParameter(rates.getParameter("alpha"));

    double lastLogLike = starLikelihood(model, rates, sequences);
    const size_t nParam = params.size();

    if(verbose) {
        std::clog << "intial: " << lastLogLike << '\n';
        std::clog.flush();
    }

    size_t iter = 0;
    for(iter = 0; iter < MAX_ROUNDS; iter++) {
        bool anyImproved = false;

        // First, substitution model
        nlopt::opt opt(nlopt::LN_BOBYQA, nParam);
        //nlopt::opt opt(nlopt::LN_COBYLA, nParam);
        NlOptParams optParams { &sequences, &params, &model, &rates };
        opt.set_max_objective(nlLogLike, &optParams);

        std::vector<double> lowerBounds(nParam, -std::numeric_limits<double>::max());
        std::vector<double> upperBounds(nParam, std::numeric_limits<double>::max());
        for(size_t i = 0; i < nParam; i++) {
            if(!params[i].hasConstraint())
                continue;

            // TODO: changes in bpp 2.1
            bpp::Interval* constraint = dynamic_cast<bpp::Interval*>(params[i].getConstraint());
            assert(constraint != nullptr);
            lowerBounds[i] = constraint->getLowerBound();
            if(lowerBounds[i] == 0.0) lowerBounds[i] += 1e-7;
            upperBounds[i] = constraint->getUpperBound();
        }
        opt.set_lower_bounds(lowerBounds);
        opt.set_upper_bounds(upperBounds);
        opt.set_ftol_abs(IMPROVE_THRESH);
        opt.set_maxeval(MAX_ITER);

        std::vector<double> x(nParam);
        for(size_t i = 0; i < nParam; i++) {
            x[i] = params[i].getValue();
        }

        double logLike;
        try {
            const int nlOptResult = opt.optimize(x, logLike);
            std::clog << "Optimization finished with " << nlOptResult << '\n';
        } catch(std::exception& e) {
            std::clog << "Optimization failed.\n";
        }
        //if(nlOptResult != nlopt::SUCCESS)

        for(size_t i = 0; i < nParam; i++) {
            params[i].setValue(x[i]);
        }

        if(verbose) {
            std::clog << "iteration " << iter << ": " << lastLogLike << " ->\t" << logLike << '\t' << logLike - lastLogLike << '\n';
            std::clog.flush();
        }

        if(std::abs(logLike - lastLogLike) > IMPROVE_THRESH)
            anyImproved = true;
        lastLogLike = logLike;

        // then, branch lengths
        estimateBranchLengths(model,
                              rates,
                              sequences);
        logLike = starLikelihood(model, rates, sequences);

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

    if(verbose) {
        for(std::size_t i = 0; i < params.size(); i++) {
            std::clog << params[i].getName() << "=" << params[i].getValue() << "\n";
        }
    }

    return iter;
}

} // namespace
