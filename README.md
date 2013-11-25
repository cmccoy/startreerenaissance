# Requirements

    libprotobuf-dev
    beagle-2.1
    bpp-core-dev 2.0
    bpp-phyl-dev 2.0
    bpp-seq-dev 2.0
    libz-dev

# Programs

This package fits the GTR model to pairwise alignments in BAM files.

Functionality is separated into two parts:

    * `build_mutation_matrices` transforms each aligned sequence in a BAM file into a 4x4 matrix containing counts of each nucleotide substitution
    * `fit_star` fits the GTR model the output of `build_mutation_matrices`, outputting a JSON document with the results.
