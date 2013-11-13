set(SAM_FILES
    bgzf.c kstring.c bam_aux.c bam.c bam_import.c sam.c bam_index.c
    bam_pileup.c bam_lpileup.c bam_md.c razf.c faidx.c bedidx.c
    knetfile.c bam_sort.c sam_header.c bam_reheader.c kprobaln.c bam_cat.c)

set(CMAKE_C_FLAGS "-g -O2 ${CMAKE_C_FLAGS}")

add_library(bam STATIC ${SAM_FILES})
set_property(TARGET bam PROPERTY COMPILE_DEFINITIONS
             _FILE_OFFSET_BITS=64 _LARGEFILE64_SOURCE _USE_KNETFILE)
target_link_libraries(bam pthread z c)
