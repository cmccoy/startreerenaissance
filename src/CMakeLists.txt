# Protocol buffers
include(FindProtobuf)
PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders mutationio.proto)

add_library(libgtrfit
            mutationio.pb.cc
            gtr.cpp)
set_target_properties(libgtrfit PROPERTIES OUTPUT_NAME gtrfit)
target_link_libraries(libgtrfit protobuf z c pthread bam)
set_property(TARGET libgtrfit PROPERTY COMPILE_FLAGS "-fopenmp")
set_property(TARGET libgtrfit PROPERTY LINK_FLAGS "-fopenmp")

# Build mutation matrices
add_executable(build_mutation_matrices
  build_mutation_matrices.cpp
  mutationio.pb.cc)
target_link_libraries(build_mutation_matrices libgtrfit)

# GTR fitting
add_executable(fit_gtr fit_gtr.cpp)
target_link_libraries(fit_gtr boost-programoptions libgtrfit)
set_property(TARGET fit_gtr PROPERTY COMPILE_FLAGS "-fopenmp")
set_property(TARGET fit_gtr PROPERTY LINK_FLAGS "-fopenmp")
