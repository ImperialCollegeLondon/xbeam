cmake_minimum_required(VERSION 3.9)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")
enable_language(Fortran)
project(xbeam_lib)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(BUILD_SHARED_LIBS ON)
file(GLOB SOURCES src/xbeam_base/*.f90 src/*.f90)
add_library(xbeam SHARED ${SOURCES})

find_package(MKL)
if (MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIR})
    target_link_libraries(xbeam PUBLIC ${MKL_LIBRARIES})
else()
    find_package(BLAS)
    target_link_libraries(xbeam PUBLIC ${BLAS_LIBRARIES})
    find_package(LAPACK REQUIRED)
    target_link_libraries(xbeam PUBLIC ${LAPACK_LIBRARIES})
endif()

FIND_PACKAGE(OpenMP)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(xbeam PUBLIC OpenMP::OpenMP_Fortran)
endif()

# Add debug flags
## For all compilers
target_compile_options(xbeam PUBLIC $<$<CONFIG:DEBUG>:>)

# Add release flags
## For all compilers
target_compile_options(xbeam PUBLIC $<$<CONFIG:RELEASE>:-funroll-loops>)

# Custom flags for:
## For Intel Fortran
## To do - test optimisation flags for IFORT before implementation
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_compile_options(xbeam PUBLIC $<$<CONFIG:DEBUG>:-warn all>)
    target_compile_options(xbeam PUBLIC $<$<CONFIG:RELEASE>:-heap-arrays>)
endif()

## For GNU gfortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(xbeam PUBLIC $<$<CONFIG:DEBUG>:-Wall>)
    target_compile_options(xbeam PUBLIC $<$<CONFIG:RELEASE>:-ofast>)
endif()

install(TARGETS xbeam
    LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)
