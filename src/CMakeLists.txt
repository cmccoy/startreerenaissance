# Protocol buffers
include(FindProtobuf)
PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders mutationio.proto)

# Build mutation matrices
add_executable(build_mutation_matrices
  build_mutation_matrices.cpp
  mutationio.pb.cc)
target_link_libraries(build_mutation_matrices protobuf z c pthread bam)

# GTR fitting
add_executable(fit_gtr
  fit_gtr.cpp
  gtr.cpp
  mutationio.pb.cc)
target_link_libraries(fit_gtr protobuf z boost-programoptions)
set_property(TARGET fit_gtr PROPERTY COMPILE_FLAGS "-fopenmp")
set_property(TARGET fit_gtr PROPERTY LINK_FLAGS "-fopenmp")
