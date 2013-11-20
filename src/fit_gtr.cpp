#include <algorithm>
#include <vector>
#include <fstream>
#include <iostream>
#include "mutationio.pb.h"

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/gzip_stream.h>

#include <boost/program_options.hpp>
#include <json/json.h>
#include <json/value.h>

#include "gtr.hpp"
#include "sequence.hpp"

using namespace gtr;
using Eigen::Matrix4d;
namespace po = boost::program_options;

std::vector<Sequence> loadSequencesFromFile(const std::string& path)
{
    std::fstream in(path, std::ios::in | std::ios::binary);
    std::clog << "Loading from " << path << '\n';
    assert(in.good() && "Input stream is not good.");
    google::protobuf::io::IstreamInputStream raw_in(&in);
    google::protobuf::io::GzipInputStream zip_in(&raw_in);

    std::vector<Sequence> sequences;

    while(true) {
        google::protobuf::io::CodedInputStream coded_in(&zip_in);
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
        Sequence sequence;
        if(m.has_name())
            sequence.name = m.name();
        for(size_t i = 0; i < 4; i++)
            for(size_t j = 0; j < 4; j++)
                sequence.substitutions(i, j) = m.mutations(4 * i + j);
        if(m.has_distance())
            sequence.distance = m.distance();
        else {
            double d = 1 - sequence.substitutions.diagonal().sum() / sequence.substitutions.sum();
            if(d == 0)
                d = 1e-3;
            sequence.distance = d;
        }
        sequences.push_back(std::move(sequence));
    }
    return sequences;
}

void writeResults(std::ostream& out,
                  const GTRParameters parameters,
                  const std::vector<Sequence>& sequences)
{
    Json::Value root;
    Json::Value modelNode(Json::objectValue);
    Json::Value thetaNode(Json::arrayValue);
    Json::Value piNode(Json::arrayValue);
    std::vector<std::string> names {"ct", "at", "gt", "ac", "cg", "ag"};
    for(size_t i = 0; i < parameters.params.size(); i++)
        modelNode[names[i]] = parameters.params[i];
    for(size_t i = 0; i < parameters.theta.size(); i++)
        thetaNode.append(parameters.theta[i]);
    modelNode["theta"] = thetaNode;

    const Eigen::Vector4d pi = thetaToPi(parameters.theta);
    for(size_t i = 0; i < pi.size(); i++)
        piNode.append(pi[i]);
    modelNode["pi"] = piNode;

    // Matrices
    Json::Value qNode(Json::arrayValue);
    const Eigen::Matrix4d q = parameters.createQMatrix();
    for(size_t i = 0; i < 4; i++) {
        for(size_t j = 0; j < 4; j++) {
            qNode.append(q(i, j));
        }
    }
    modelNode["Q"] = qNode;

    auto f = [](double acc, const Sequence & s) { return acc + s.distance; };
    const double meanBranchLength = std::accumulate(sequences.begin(), sequences.end(), 0.0, f) / sequences.size();
    const Eigen::Matrix4d p = parameters.createModel().createPMatrix(meanBranchLength);
    root["mean_branch_length"] = meanBranchLength;

    Json::Value pNode(Json::arrayValue);
    for(size_t i = 0; i < 4; i++) {
        for(size_t j = 0; j < 4; j++) {
            pNode.append(p(i, j));
        }
    }
    modelNode["P_mean"] = pNode;

    root["model"] = modelNode;

    root["logLikelihood"] = starLikelihood(parameters.createModel(), sequences);
    Json::Value blNode(Json::arrayValue);
    for(const Sequence& sequence : sequences)
        blNode.append(sequence.distance);
    root["branch_lengths"] = blNode;

    out << root << '\n';
}

int main(const int argc, const char** argv)
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    // command-line parsing
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "Produce help message")
    ("input-file,i", po::value<std::string>(), "input file [required]")
    ("output-file,o", po::value<std::string>(), "output file [required]");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(desc).run(), vm);
    po::notify(vm);

    if(vm.count("help") || vm.count("h")) {
        std::clog << desc << '\n';
        return 1;
    }

    if(!vm.count("input-file")) {
        std::clog << "missing input.\n";
        return 1;
    }

    if(!vm.count("output-file")) {
        std::clog << "missing output.\n";
        return 1;
    }

    std::vector<Sequence> sequences = loadSequencesFromFile(vm["input-file"].as<std::string>());

    std::cout << sequences.size() << " sequences." << '\n';

    Matrix4d m = std::accumulate(sequences.begin() + 1, sequences.end(), sequences[0].substitutions,
                                 [](const Matrix4d& acc, const Sequence& s) { return acc + s.substitutions; });
    std::cout << "Substitution counts[raw]:\n " << m << '\n';


    gtr::GTRParameters params;
    empiricalModel(sequences, params);

    gtr::GTRModel model = params.createModel();

    std::cout << "Initial log-like: " << starLikelihood(model, sequences) << '\n';

    optimize(params, sequences);

    auto f = [](double acc, const Sequence & s) { return acc + s.distance; };
    const double meanBranchLength = std::accumulate(sequences.begin(), sequences.end(), 0.0, f) / sequences.size();

    std::cout << "Mean branch length: " << meanBranchLength << '\n';

    std::cout << "Final log-like: " << starLikelihood(params.createModel(), sequences) << '\n';

    std::cout << "ct at gt ac cg ag = " << params.params.transpose() << '\t' << 1 << '\n';
    std::cout << "pi = " << params.createBaseFrequencies() << '\n';
    std::cout << "Q=\n" << params.createQMatrix() << '\n';
    std::cout << "P(" << meanBranchLength << ")=\n" << params.createModel().createPMatrix(meanBranchLength) << '\n';

    std::ofstream out(vm["output-file"].as<std::string>());
    writeResults(out, params, sequences);

    google::protobuf::ShutdownProtobufLibrary();
    return 0;
}
