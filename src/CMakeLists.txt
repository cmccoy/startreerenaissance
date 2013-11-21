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

add_library(libgtrfit
            mutationio.pb.cc
            star_optim.cpp
            gtr.cpp)
set_target_properties(libgtrfit PROPERTIES OUTPUT_NAME gtrfit)
target_link_libraries(libgtrfit protobuf z c pthread bam nlopt bpp-phyl bpp-core bpp-seq
                      ${HMS_BEAGLE_LIBRARIES})
#set_property(TARGET libgtrfit PROPERTY COMPILE_FLAGS "-fopenmp")
#set_property(TARGET libgtrfit PROPERTY LINK_FLAGS "-fopenmp")

# Build mutation matrices
add_executable(build_mutation_matrices
  build_mutation_matrices.cpp
  mutationio.pb.cc)
target_link_libraries(build_mutation_matrices libgtrfit)

# GTR fitting
add_executable(fit_gtr fit_gtr.cpp)
target_link_libraries(fit_gtr boost-programoptions libgtrfit)
target_link_libraries(fit_gtr json)
target_link_libraries(fit_gtr nlopt)
set_property(TARGET fit_gtr PROPERTY COMPILE_FLAGS "-fopenmp")
set_property(TARGET fit_gtr PROPERTY LINK_FLAGS "-fopenmp")
