cmake_minimum_required(VERSION 3.8)
project(eScatter CXX CUDA)

option(WITH_UTILS "Build all utilities." OFF)
option(WITH_DOCS "Build documentation." OFF)
option(WITH_TEST "Build all tests." OFF)
option(WITH_CSTOOL "Build cstool (deprecated)." OFF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 11)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

string(APPEND CMAKE_CXX_FLAGS " -O3")
#string(APPEND CMAKE_CXX_FLAGS " -Wall -ggdb") #debug
string(APPEND CMAKE_CUDA_FLAGS " -O3 -prec-sqrt=false")
#string(APPEND CMAKE_CUDA_FLAGS " -Xptxas=-v --compiler-options -Wall") #debug
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_37,code=compute_37")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_52,code=compute_52")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_60,code=compute_60")
string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_61,code=compute_61") #TUD HPC25
if(CMAKE_CUDA_COMPILER_VERSION LESS 9.0)
	string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_20,code=compute_20") #TUD T7500
endif(CMAKE_CUDA_COMPILER_VERSION LESS 9.0)

if(WITH_DOCS)
	find_package(Doxygen REQUIRED)
	add_custom_target(doc ALL ${DOXYGEN_EXECUTABLE} ${PROJECT_SOURCE_DIR}/Doxyfile COMMENT "Generating documentation" VERBATIM)
endif(WITH_DOCS)

if(WITH_TEST)
	enable_testing()
	find_package(GTest REQUIRED)
	add_subdirectory(test)
endif(WITH_TEST)

include_directories(src 3rdparty/cub 3rdparty/GSL)
add_subdirectory(src)
