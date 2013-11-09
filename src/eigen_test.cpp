#include <atomic>
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
            return - model.logLikelihood(s);
        };
        boost::uintmax_t max_iter = 50;
        s.distance = brent_find_minima(f, 1e-7, 1.0, 8, max_iter).second;
    }
}

Eigen::Vector4d count_base_frequences(const std::vector<gtr::Sequence>& sequences)
{
    Eigen::Vector4d result;

    for(const gtr::Sequence& s : sequences)
        result += s.transitions.colwise().sum();

    result /= result.sum();
    return result;
}

int main()
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;
    std::vector<gtr::Sequence> sequences;

    std::fstream in("test.bin.gz", std::ios::in | std::ios::binary);
    google::protobuf::io::IstreamInputStream raw_in(&in);
    google::protobuf::io::GzipInputStream zip_in(&raw_in);
    google::protobuf::io::CodedInputStream coded_in(&zip_in);

    while(in.good()) {
        uint32_t size;
        bool success = false;
        success = coded_in.ReadVarint32(&size);
        if(!success) break;
        mutationio::MutationCount m;
        m.ParseFromBoundedZeroCopyStream(&raw_in, size);

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

    return 0;
}
