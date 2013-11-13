cmake_minimum_required (VERSION 2.8.3)
project (fit_gtr)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/lib/cmake/cmake-modules")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/lib/cmake/cxx11/Modules")

include(GetGitRevisionDescription)
git_describe(GIT_DESCRIPTION --always --long)
configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_SOURCE_DIR}/src/config.h")

include(FindCXXFeatures)

if(NOT _HAS_CXX11_FLAG)
   message( FATAL_ERROR "-std=c++11 flag not supported." )
endif()

include_directories("src")
include_directories("lib/eigen_3.2.0")
include_directories("lib/samtools")
include_directories("lib/boost_1.54.0")
set(CMAKE_CXX_FLAGS "-g -Wall -Wextra -O2 ${CXX11_COMPILER_FLAGS}")

add_subdirectory(lib/samtools)
add_subdirectory(lib/boost_1.54.0/libs/program_options)
add_subdirectory(src)
