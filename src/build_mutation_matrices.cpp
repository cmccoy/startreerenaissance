#include "sam.h"
#include "faidx.h"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/gzip_stream.h>

#include "mutationio.pb.h"

namespace po = boost::program_options;
namespace protoio = google::protobuf::io;

bool endsWith(const std::string& s, const std::string& suffix)
{
    if(s.length() >= suffix.length()) {
        return (0 == s.compare(s.length() - suffix.length(), suffix.length(), suffix));
    } else {
        return false;
    }
}

inline int nt16ToIdx(const int b)
{
    switch(b) {
        case 1: return 0; // A
        case 2: return 1; // C
        case 4: return 2; // G
        case 8: return 3; // T
        default: return 4; // N
    }
}

int usage(po::options_description& desc)
{
    std::cerr << "Usage: build_mutation_matrices [options] <ref.fasta> <in.bam> <out.bin>\n";
    std::cerr << desc << '\n';
    return 1;
}


struct SamFile {
    SamFile(const std::string& path, const std::string& mode = "rb", void* extra = nullptr) :
        fp(samopen(path.c_str(), mode.c_str(), extra))
    {
        assert(fp != nullptr && "Failed to open BAM");
    }

    ~SamFile()
    {
        if(fp != nullptr)
            samclose(fp);
    }

    samfile_t* fp;
};

int main(int argc, char* argv[])
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;

    std::string fastaPath, bamPath, outputPath;

    bool no_ambiguous = false;
    size_t maxRecords = 0;

    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h", "Produce help message")
    ("no-ambiguous", po::bool_switch(&no_ambiguous), "Do not include ambiguous sites")
    ("max-records,n", po::value<size_t>(&maxRecords), "Maximum number of records to parse")
    ("input-fasta,f", po::value<std::string>(&fastaPath)->required(), "Path to (indexed) FASTA file")
    ("input-bam,i", po::value<std::string>(&bamPath)->required(), "Path to BAM")
    ("output-file,o", po::value<std::string>(&outputPath)->required(), "Path to output file");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help")) {
        std::cout << desc << '\n';
        return 0;
    }

    po::notify(vm);

    faidx_t* fidx = fai_load(fastaPath.c_str());
    assert(fidx != NULL && "Failed to load FASTA index");

    SamFile in(bamPath);

    bam1_t* b = bam_init1();
    size_t processed = 0;

    std::vector<std::string> targetBases(in.fp->header->n_targets);
    std::vector<int> targetLen(in.fp->header->n_targets);

    std::fstream out(outputPath, std::ios::out | std::ios::trunc | std::ios::binary);
    protoio::OstreamOutputStream rawOut(&out);
    protoio::GzipOutputStream zipOut(&rawOut);
    protoio::ZeroCopyOutputStream* outptr = &rawOut;
    if(endsWith(outputPath, ".gz"))
        outptr = &zipOut;
    protoio::CodedOutputStream codedOut(outptr);

    while(samread(in.fp, b) >= 0) {
        if(maxRecords > 0 && processed > maxRecords)
            break;
        processed++;
        const char* target_name = in.fp->header->target_name[b->core.tid];
        if(targetBases[b->core.tid].empty()) {
            char* ref = fai_fetch(fidx, target_name, &targetLen[b->core.tid]);
            assert(ref != nullptr && "Missing reference");
            targetBases[b->core.tid] = ref;
            free(ref);
        }

        const uint32_t* cigar = bam1_cigar(b);
        const uint8_t* seq = bam1_seq(b);
        uint32_t qi = 0, ri = b->core.pos;
        const std::string& ref = targetBases[b->core.tid];

        std::vector<int> mutations(16, 0);

        const int8_t* bq = reinterpret_cast<int8_t*>(bam_aux_get(b, "bq"));
        if(no_ambiguous) {
            assert(bq != NULL && "No bq tag");
        }
        for(uint32_t cidx = 0; cidx < b->core.n_cigar; cidx++) {
            const uint32_t clen = bam_cigar_oplen(cigar[cidx]);
            const uint32_t consumes = bam_cigar_type(cigar[cidx]); // bit 1: consume query; bit 2: consume reference
            if((consumes & 0x3) == 0x3) {  // Reference and query
                for(uint32_t i = 0; i < clen; i++) {
                    const int qb = nt16ToIdx(bam1_seqi(seq, qi + i)),
                              rb = nt16ToIdx(bam_nt16_table[static_cast<int>(ref[ri + i])]);
                    if(qb < 4 && rb < 4 && (!no_ambiguous || bq[qi + i] % 100 == 0)) {
                        mutations[(rb * 4) + qb] += 1;
                    }
                }
            }
            if(consumes & 0x1) // Consumes query
                qi += clen;
            else if(consumes & 0x2) // Consumes reference
                ri += clen;
        }

        mutationio::MutationCount count;
        count.set_name(bam1_qname(b));
        for(const int i : mutations)
            count.add_mutations(i);

        codedOut.WriteVarint32(count.ByteSize());
        count.SerializeWithCachedSizes(&codedOut);
    }

    bam_destroy1(b);
    fai_destroy(fidx);

    google::protobuf::ShutdownProtobufLibrary();

    return 0;
}
