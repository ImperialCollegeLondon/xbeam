cmake_minimum_required(VERSION 3.5)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
enable_language(Fortran)
project(xbeam_lib)

set(CMAKE_VERBOSE_MAKEFILE ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_F90_COMPILER gfortran)

file(GLOB SOURCES src/xbeam_base/*.f90 src/*.f90)
add_library(xbeam SHARED ${SOURCES})



set(BUILD_SHARED_LIBS ON)
find_package(MKL)
if (MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIR})
    target_link_libraries(xbeam PUBLIC ${MKL_LIBRARIES})
else()
    set(BLA_STATIC ON)
    find_package(LAPACK)
    if(NOT LAPACK_FOUND)
        set(BLA_STATIC ON)
        set(BLA_VENDOR Generic)
        find_package(LAPACK REQUIRED)
    endif()
    target_link_libraries(xbeam PUBLIC ${LAPACK_LIBRARIES})
endif()

FIND_PACKAGE(OpenMP)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(xbeam PUBLIC OpenMP::OpenMP_Fortran)
endif()


# Add debug flags
## For all compilers
target_compile_options(xbeam PUBLIC $<$<CONFIG:DEBUG>: >)

# Add release flags
## For all compilers
target_compile_options(xbeam PUBLIC $<$<CONFIG:RELEASE>: -funroll-loops>)

# Custom flags for:
## For Intel Fortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    target_compile_options(xbeam PUBLIC $<$<CONFIG:DEBUG>:-warn all>)
    target_compile_options(xbeam PUBLIC $<$<CONFIG:RELEASE>:-heap-arrays>)
endif()

## For GNU gfortran
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(xbeam PUBLIC $<$<CONFIG:DEBUG>:-Wall>)
    target_compile_options(xbeam PUBLIC $<$<CONFIG:RELEASE>:-ftree-parallelize-loops=4>)
endif()


install(TARGETS xbeam
    LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/lib)
