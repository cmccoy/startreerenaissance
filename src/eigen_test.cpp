#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include <boost/math/tools/minima.hpp>
#include "mutationlist.pb.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/gzip_stream.h>

#include "gtr.hpp"

using namespace gtr;
using boost::math::tools::brent_find_minima;

double star_likelihood(const GTRModel& model,
                       const std::vector<Sequence>& sequences)
{
    double result = 0.0;

#pragma omp parallel for reduction(+:result)
    for(size_t i = 0; i < sequences.size(); i++) {
        result += model.logLikelihood(sequences[i]);
    }
    return result;
}

void estimate_branch_lengths(const GTRModel& model,
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
        std::pair<double, double> res =  brent_find_minima(f, 1e-9, 1.0, 50, max_iter);
        s.distance = res.first;
    }
}

Eigen::Vector4d count_base_frequences(const std::vector<gtr::Sequence>& sequences)
{
    Eigen::Vector4d result;
    result.fill(1);

    for(const gtr::Sequence& s : sequences)
        result += s.transitions.colwise().sum();

    result /= result.sum();
    return result;
}

void empirical_model(const std::vector<gtr::Sequence>& sequences,
                     gtr::GTRParameters& model)
{
    model.pi = count_base_frequences(sequences);

    Matrix4d result;
    result.fill(0);

    for(const gtr::Sequence& s : sequences) {
        result += s.transitions;
    }

    model.params << result(0, 1), result(0, 2), result(0, 3),
                    result(1, 2), result(1, 3), result(2, 3);
    model.params /= model.params[5];

    std::cout << model.params.transpose() << '\n'
              << model.pi.transpose()  << '\n' << result << '\n';
}

std::vector<gtr::Sequence> load_sequences_from_file(const std::string& path)
{
    std::fstream in(path, std::ios::in | std::ios::binary);
    google::protobuf::io::IstreamInputStream raw_in(&in);
    google::protobuf::io::GzipInputStream zip_in(&raw_in);
    google::protobuf::io::CodedInputStream coded_in(&zip_in);

    std::vector<gtr::Sequence> sequences;

    while(true) {
        uint32_t size = 0;
        bool success = false;
        success = coded_in.ReadVarint32(&size);
        if(!success) break;
        mutationio::MutationCount m;
        std::string s;
        coded_in.ReadString(&s, size);
        success = m.ParseFromString(s);
        assert(success && "Failed to parse");
        assert(m.mutations_size() == 16 && "Unexpected mutations count");
        gtr::Sequence sequence;
        sequence.distance = m.has_distance() ? m.distance() : 0.1;
        if(m.has_name())
            sequence.name = m.name();
        for(size_t i = 0; i < 4; i++)
            for(size_t j = 0; j < 4; j++)
                sequence.transitions(i, j) = m.mutations(4*i + j);
        sequences.push_back(std::move(sequence));
    }
    return sequences;
}

int main()
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;
    std::vector<gtr::Sequence> sequences = load_sequences_from_file("test.muts.pb.gz");

    std::cout << sequences.size() << " sequences." << '\n';

    google::protobuf::ShutdownProtobufLibrary();

    gtr::GTRParameters params;
    empirical_model(sequences, params);

    gtr::GTRModel model = params.buildModel();

    std::cout << "Initial log-like: " << star_likelihood(model, sequences) << '\n';

    estimate_branch_lengths(model, sequences);

    auto f = [](double acc, const Sequence& s) { return acc + s.distance; };
    const double mean_branch_length = std::accumulate(sequences.begin(), sequences.end(), 0.0, f) / sequences.size();
    std::cout << "Mean branch length: " << mean_branch_length << '\n';
    std::cout << "Final log-like: " << star_likelihood(model, sequences) << '\n';

    return 0;
}
