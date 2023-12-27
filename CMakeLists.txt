### To do:
# 1. On M1-Macs, do:
# cmake -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++ -DCMAKE_C_COMPILER=$(brew --prefix llvm)/bin/clang -DCVQEC_USE_OPENMP=OFF ../
###

cmake_minimum_required(VERSION 3.14)

project(qecc VERSION 0.1.0)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

################################################################################
#                                  Set options
################################################################################

option(QECC_ENABLE_WARNINGS "Enable warnings" ON)
option(QECC_USE_OPENMP "Use OpenMP" OFF)

add_library(qecc_compile_options INTERFACE)

if(QECC_ENABLE_WARNINGS)
    target_compile_options(qecc_compile_options INTERFACE -Wall;-Wextra;-Werror)
endif()


################################################################################
#                                Set depdencies
################################################################################

add_library(qecc_external_libs INTERFACE)

find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

target_link_libraries(qecc_external_libs INTERFACE ${lapackblas_libraries})

if (QECC_USE_OPENMP)
    message(STATUS "Use OpenMP")
    find_package(OpenMP REQUIRED)
    target_link_libraries(qecc_external_libs INTERFACE OpenMP::OpenMP_C)
endif()

################################################################################
#                                Set build
################################################################################

set(SOURCE_FILES
    src/checks.c
    src/printfuns.c
    src/qecc.c
    src/utils.c
    src/pmets.c
    src/sampling.c
    src/decode.c
    src/rand.c
    src/mt19937/mt19937ar.c
    src/reps.c
    src/benchmark.c
    src/hybrid.c
    src/constants.c
    src/zeigh.c
    src/logmetrics.c
    src/effective.c
    src/memory.c
    src/main.c
    src/linalg.c
)

add_library(qecc-lib STATIC
    ${SOURCE_FILES}
)

target_include_directories(cvqec-lib PUBLIC "${PROJECT_SOURCE_DIR}/include")

add_executable(qecc
    src/main.c
)
target_link_libraries(qecc PUBLIC qecc-lib qecc_compile_options)