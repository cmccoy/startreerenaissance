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

#include "star_optim.hpp"
#include "sequence.hpp"

// Bio++
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Phyl/Model/GTR.h>
#include <Bpp/Seq/Alphabet/DNA.h>

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

        // TODO: use mutation count
        //double d = 1 - sequence.substitutions.diagonal().sum() / sequence.substitutions.sum();
        //if(d == 0)
        //d = 1e-3;
        sequence.distance = 0.1;

        sequences.push_back(std::move(sequence));
    }
    return sequences;
}

void writeResults(std::ostream& out,
                  const bpp::SubstitutionModel& model,
                  const bpp::DiscreteDistribution& rates,
                  const std::vector<Sequence>& sequences,
                  const bool include_branch_lengths = true)
{
    Json::Value root;
    Json::Value modelNode(Json::objectValue);
    Json::Value rateNode(Json::objectValue);
    Json::Value parameterNode(Json::objectValue);
    Json::Value piNode(Json::arrayValue);
    std::vector<std::string> names {"a", "b", "c", "d", "e", "ag"};
    bpp::ParameterList p = model.getParameters();
    for(size_t i = 0; i < p.size(); i++) {
        parameterNode[p[i].getName()] = p[i].getValue();
    }
    modelNode["name"] = model.getName();
    modelNode["parameters"] = parameterNode;

    //rateNode["name"] = rates.getName();
    p = rates.getParameters();
    for(size_t i = 0; i < p.size(); i++) {
        rateNode[p[i].getName()] = p[i].getValue();
    }
    root["rate"] = rateNode;


    for(size_t i = 0; i < model.getNumberOfStates(); i++) {
        piNode.append(model.freq(i));
    }
    modelNode["pi"] = piNode;

    auto f = [](double acc, const Sequence & s) { return acc + s.distance; };
    const double meanBranchLength = std::accumulate(sequences.begin(), sequences.end(), 0.0, f) / sequences.size();
    root["mean_branch_length"] = meanBranchLength;

    // Matrices
    Json::Value qNode(Json::arrayValue);
    //Json::Value sNode(Json::arrayValue);
    Json::Value pNode(Json::arrayValue);

    const auto& pt = model.getPij_t(meanBranchLength);
    for(size_t i = 0; i < model.getNumberOfStates(); i++) {
        for(size_t j = 0; j < model.getNumberOfStates(); j++) {
            qNode.append(model.Qij(i, j));
            //sNode.append(model.Sij(i, j));
            pNode.append(pt(static_cast<unsigned int>(i), static_cast<unsigned int>(j)));
        }
    }
    modelNode["Q"] = qNode;
    //modelNode["S"] = sNode;
    modelNode["P_mean"] = pNode;

    root["model"] = modelNode;

    root["logLikelihood"] = star_optim::starLikelihood(model, rates, sequences);
    if(include_branch_lengths) {
        Json::Value blNode(Json::arrayValue);
        for(const Sequence& sequence : sequences)
            blNode.append(sequence.distance);
        root["branch_lengths"] = blNode;
    }

    out << root << '\n';
}

int main(const int argc, const char** argv)
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    std::string input_path, output_path;
    bool no_branch_lengths = false;

    // command-line parsing
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "Produce help message")
    ("input-file,i", po::value(&input_path), "input file [required]")
    ("output-file,o", po::value(&output_path), "output file [required]")
    ("no-branch-lengths", po::bool_switch(&no_branch_lengths), "*do not* include fit branch lengths in output");

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

    auto sumLength = [](const Eigen::Matrix4d & acc, const Sequence & s) {
        return acc + s.substitutions;
    };
    Eigen::Matrix4d m = std::accumulate(sequences.begin() + 1, sequences.end(), sequences[0].substitutions,
                                        sumLength);
    std::cout << "Substitution counts[raw]:\n " << m << '\n';

    bpp::DNA dna;
    bpp::GTR gtr(&dna);
    bpp::GammaDiscreteDistribution rates(4);
    //bpp::ConstantDistribution rates(1.0);

    std::cout << "Initial log-like: " << star_optim::starLikelihood(gtr, rates, sequences) << '\n';

    star_optim::optimize(gtr, rates, sequences);

    std::cout << "final log-like: " << star_optim::starLikelihood(gtr, rates, sequences) << '\n';

    std::ofstream out(vm["output-file"].as<std::string>());
    writeResults(out, gtr, rates, sequences, !no_branch_lengths);

    google::protobuf::ShutdownProtobufLibrary();
    return 0;
}
