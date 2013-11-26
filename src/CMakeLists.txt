# Protocol buffers
include(FindProtobuf)
PROTOBUF_GENERATE_CPP(ProtoSources ProtoHeaders mutationio.proto)

# BEAGLE
find_package(PkgConfig REQUIRED)
pkg_check_modules(HMS_BEAGLE hmsbeagle-2 REQUIRED)
include_directories(${HMS_BEAGLE_INCLUDE_DIRS})
link_directories(${HMS_BEAGLE_LIBRARY_DIRS})

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  add_definitions("-DHAVE_OMP")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_library(libstaropt
            mutationio.pb.cc
            star_optim.cpp)
set_target_properties(libstaropt PROPERTIES OUTPUT_NAME staropt)
target_link_libraries(libstaropt protobuf z c pthread bam nlopt bpp-phyl bpp-core bpp-seq
                      ${HMS_BEAGLE_LIBRARIES})

# Build mutation matrices
add_executable(build_mutation_matrices
  build_mutation_matrices.cpp
  mutationio.pb.cc)
target_link_libraries(build_mutation_matrices libstaropt)
target_link_libraries(build_mutation_matrices boost-programoptions)

# GTR fitting
add_executable(fit_star fit_star.cpp)
target_link_libraries(fit_star boost-programoptions libstaropt)
target_link_libraries(fit_star json)
target_link_libraries(fit_star nlopt)
