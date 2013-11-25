#include "sam.h"
#include "faidx.h"

#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/gzip_stream.h>

#include "mutationio.pb.h"

int usage()
{
    fprintf(stderr, "Usage: build_mutation_matrices [options] <ref.fasta> <in.bam> <out.bin>\n\n");
    fprintf(stderr, "Options: -n INT    Maximum number number of records\n");
    fprintf(stderr, "Options: -m        Only include sites which map unambiguously\n");
    return 1;
}

inline int nt16_to_idx(const int b)
{
    switch(b) {
        case 1: return 0; // A
        case 2: return 1; // C
        case 4: return 2; // G
        case 8: return 3; // T
        default: return 4; // N

    }
}

int main(int argc, char* argv[])
{
    GOOGLE_PROTOBUF_VERIFY_VERSION;
    if(argc < 4) {
        return usage();
    }

    size_t n = 0;
    int ambiguous = 1;
    char c;
    while((c = getopt(argc, argv, "mn:h?")) >= 0) {
        switch(c) {
            case 'n': n = atoi(optarg); break;
            case 'm': ambiguous = 0; break;
            default: return usage();
        }
    }

    argv = argv + optind;
    argc = argc - optind;

    if(argc != 3)
        return usage();

    faidx_t* fidx = fai_load(argv[0]);
    assert(fidx != NULL && "Failed to load FASTA index");

    samfile_t* in = NULL;
    if((in = samopen(argv[1], "rb", NULL)) == 0) {
        fprintf(stderr, "Failed reading %s", argv[1]);
        return 1;
    }

    bam1_t* b = bam_init1();
    size_t processed = 0;

    std::vector<std::string> target_bases(in->header->n_targets);
    std::vector<int> target_len(in->header->n_targets);

    std::fstream out(argv[2], std::ios::out | std::ios::trunc | std::ios::binary);
    google::protobuf::io::OstreamOutputStream raw_out(&out);
    google::protobuf::io::GzipOutputStream zip_out(&raw_out);
    google::protobuf::io::CodedOutputStream coded_out(&zip_out);

    while(samread(in, b) >= 0) {
        if(n > 0 && processed > n)
            break;
        processed++;
        const char* target_name = in->header->target_name[b->core.tid];
        if(target_bases[b->core.tid].empty()) {
            char* ref = fai_fetch(fidx, target_name, &target_len[b->core.tid]);
            assert(ref != nullptr && "Missing reference");
            target_bases[b->core.tid] = ref;
            free(ref);
        }

        const uint32_t* cigar = bam1_cigar(b);
        const uint8_t* seq = bam1_seq(b);
        uint32_t qi = 0, ri = b->core.pos;
        const std::string& ref = target_bases[b->core.tid];

        std::vector<int> mutations(16, 0);

        const int8_t* bq = reinterpret_cast<int8_t*>(bam_aux_get(b, "bq"));
        if(!ambiguous) {
            assert(bq != NULL && "No bq tag");
        }
        for(uint32_t cidx = 0; cidx < b->core.n_cigar; cidx++) {
            const uint32_t clen = bam_cigar_oplen(cigar[cidx]);
            const uint32_t consumes = bam_cigar_type(cigar[cidx]); // bit 1: consume query; bit 2: consume reference
            if((consumes & 0x3) == 0x3) {  // Reference and query
                for(uint32_t i = 0; i < clen; i++) {
                    const int qb = nt16_to_idx(bam1_seqi(seq, qi + i)),
                              rb = nt16_to_idx(bam_nt16_table[static_cast<int>(ref[ri + i])]);
                    if(qb < 4 && rb < 4 && (ambiguous || bq[qi + i] % 100 == 0)) {
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

        coded_out.WriteVarint32(count.ByteSize());
        count.SerializeWithCachedSizes(&coded_out);
    }

    bam_destroy1(b);
    fai_destroy(fidx);

    samclose(in);
    google::protobuf::ShutdownProtobufLibrary();

    return 0;
}
