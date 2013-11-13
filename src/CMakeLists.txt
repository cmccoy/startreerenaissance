# Protocol buffers
include(FindProtobuf)
PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders mutationio.proto)

add_library(gtrfitlib
            mutationio.pb.cc
            gtr.cpp)
target_link_libraries(gtrfitlib protobuf z c pthread bam)
set_property(TARGET gtrfitlib PROPERTY COMPILE_FLAGS "-fopenmp")
set_property(TARGET gtrfitlib PROPERTY LINK_FLAGS "-fopenmp")

# Build mutation matrices
add_executable(build_mutation_matrices
  build_mutation_matrices.cpp
  mutationio.pb.cc)
target_link_libraries(build_mutation_matrices gtrfitlib)

# GTR fitting
add_executable(fit_gtr fit_gtr.cpp)
target_link_libraries(fit_gtr boost-programoptions gtrfitlib)
set_property(TARGET fit_gtr PROPERTY COMPILE_FLAGS "-fopenmp")
set_property(TARGET fit_gtr PROPERTY LINK_FLAGS "-fopenmp")
