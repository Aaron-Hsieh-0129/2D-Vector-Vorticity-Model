cmake_minimum_required(VERSION 3.10)

# Set policy CMP0135 to NEW to avoid timestamp warnings
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Define project
project(VVM2D VERSION 0.1.0 LANGUAGES CXX Fortran)

# Allow user to specify GCC installation path
set(GCC_HOME "" CACHE STRING "Path to GCC installation")
if(GCC_HOME STREQUAL "")
    set(GCC_HOME "/home/Aaron/gcc13" CACHE STRING "Path to GCC installation" FORCE)
endif()
set(GCC_LIB_DIR "${GCC_HOME}/lib64")
set(GCC_LIB32_DIR "${GCC_HOME}/lib")

# Set GCC compilers
set(CMAKE_C_COMPILER "${GCC_HOME}/bin/gcc" CACHE STRING "Set C compiler")
set(CMAKE_CXX_COMPILER "${GCC_HOME}/bin/g++" CACHE STRING "Set C++ compiler")
set(CMAKE_Fortran_COMPILER "${GCC_HOME}/bin/gfortran" CACHE STRING "Set Fortran compiler")

# Find OpenMP
find_package(OpenMP REQUIRED)

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set where dependencies will be downloaded and installed
set(DEP_INSTALL_DIR "${CMAKE_SOURCE_DIR}/_deps/install")
file(MAKE_DIRECTORY ${DEP_INSTALL_DIR})

# Include ExternalProject for dependency management
include(ExternalProject)

# Get number of available cores
include(ProcessorCount)
ProcessorCount(N_CORES)
if(N_CORES GREATER 1)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N_CORES})
    message(STATUS "Using ${N_CORES} cores for building")
else()
    message(STATUS "Using single core for building")
endif()

# Limit cores for Boost build to avoid issues
set(BOOST_BUILD_CORES 16)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -Wall -Wno-unused-function")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -fdefault-real-8 -ffree-line-length-0")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Add OpenMP flags
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
if(OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# Set RPATH usage
set(CMAKE_INSTALL_RPATH "${DEP_INSTALL_DIR}/lib;${GCC_LIB_DIR};${GCC_LIB32_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_definitions(-DRTE_USE_CBOOL)

# Include directories
include_directories(
    include
    ${DEP_INSTALL_DIR}/include
    ${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include
    ${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include_test
)

# ExternalProject for HDF5
ExternalProject_Add(
    hdf5_ext
    URL https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.gz
    PREFIX ${CMAKE_SOURCE_DIR}/_deps/hdf5
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${DEP_INSTALL_DIR} --enable-cxx --enable-fortran --enable-shared --with-pic CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${CMAKE_Fortran_COMPILER}
    BUILD_COMMAND make -j${N_CORES}
    INSTALL_COMMAND make install
    UPDATE_COMMAND ""
    BUILD_ALWAYS FALSE
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# ExternalProject for NetCDF-C (depends on HDF5)
ExternalProject_Add(
    netcdf-c_ext
    URL https://github.com/Unidata/netcdf-c/archive/v4.3.3.1.tar.gz
    PREFIX ${CMAKE_SOURCE_DIR}/_deps/netcdf-c
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${DEP_INSTALL_DIR} --enable-netcdf-4 --enable-shared --with-pic --disable-dap CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${CMAKE_Fortran_COMPILER} CPPFLAGS=-I${DEP_INSTALL_DIR}/include LDFLAGS=-L${DEP_INSTALL_DIR}/lib
    BUILD_COMMAND make -j${N_CORES}
    INSTALL_COMMAND make install
    UPDATE_COMMAND ""
    BUILD_ALWAYS FALSE
    DEPENDS hdf5_ext
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# ExternalProject for NetCDF-CXX (depends on NetCDF-C)
ExternalProject_Add(
    netcdf-cxx4_ext
    URL https://github.com/Unidata/netcdf-cxx4/archive/v4.3.0.tar.gz
    PREFIX ${CMAKE_SOURCE_DIR}/_deps/netcdf-cxx4
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${DEP_INSTALL_DIR} --enable-shared --with-pic CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CPPFLAGS=-I${DEP_INSTALL_DIR}/include LDFLAGS=-L${DEP_INSTALL_DIR}/lib
    BUILD_COMMAND make -j${N_CORES}
    INSTALL_COMMAND make install
    UPDATE_COMMAND ""
    BUILD_ALWAYS FALSE
    DEPENDS netcdf-c_ext
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# ExternalProject for Boost
ExternalProject_Add(
    boost_ext
    URL https://archives.boost.io/release/1.87.0/source/boost_1_87_0.tar.gz
    PREFIX ${CMAKE_SOURCE_DIR}/_deps/boost
    CONFIGURE_COMMAND <SOURCE_DIR>/bootstrap.sh --prefix=${DEP_INSTALL_DIR}
    BUILD_COMMAND <SOURCE_DIR>/b2 -j${BOOST_BUILD_CORES} install --prefix=${DEP_INSTALL_DIR}
    BUILD_IN_SOURCE TRUE
    INSTALL_COMMAND ""
    UPDATE_COMMAND ""
    BUILD_ALWAYS FALSE
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)

# Add rte-rrtmgp-cpp as a subdirectory
add_subdirectory(external/rte-rrtmgp-cpp)

# Force rte-rrtmgp to use only _deps/install headers
set_target_properties(rte_rrtmgp PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${DEP_INSTALL_DIR}/include;${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include"
)

# Create the P3 Fortran shared library
add_library(p3 SHARED external/P3-microphysics/microphy_p3.f90)
set_target_properties(p3 PROPERTIES
    OUTPUT_NAME "p3"
    Fortran_PREPROCESS ON
)
if(OpenMP_Fortran_FOUND)
    target_link_libraries(p3 PRIVATE OpenMP::OpenMP_Fortran)
endif()

# Source files for VVM2D
file(GLOB_RECURSE SRC src/*.cpp src/*.hpp)
add_executable(vvm2d ${SRC})

# Link Fortran library and runtime to the executable
target_link_libraries(vvm2d PRIVATE p3 gfortran)

# Link OpenMP libraries
if(OpenMP_CXX_FOUND)
    target_link_libraries(vvm2d PRIVATE OpenMP::OpenMP_CXX)
endif()

# Specify expected library paths (will be built by ExternalProject)
set(HDF5_LIB ${DEP_INSTALL_DIR}/lib/libhdf5.so)
set(HDF5_CXX_LIB ${DEP_INSTALL_DIR}/lib/libhdf5_cpp.so)
set(NETCDF_C_LIB ${DEP_INSTALL_DIR}/lib/libnetcdf.so)
set(NETCDF_CXX_LIB ${DEP_INSTALL_DIR}/lib/libnetcdf_c++4.so)
set(BOOST_SYSTEM_LIB ${DEP_INSTALL_DIR}/lib/libboost_system.so)
set(BOOST_FILESYSTEM_LIB ${DEP_INSTALL_DIR}/lib/libboost_filesystem.so)

# Link libraries
target_link_libraries(vvm2d PRIVATE
    rte_rrtmgp
    ${HDF5_LIB}
    ${HDF5_CXX_LIB}
    ${NETCDF_C_LIB}
    ${NETCDF_CXX_LIB}
    ${BOOST_SYSTEM_LIB}
    ${BOOST_FILESYSTEM_LIB}
    stdc++
)

# Add dependencies
add_dependencies(vvm2d hdf5_ext netcdf-c_ext netcdf-cxx4_ext boost_ext)

# Set runtime library paths (-rpath)
set_target_properties(vvm2d PROPERTIES
    LINK_FLAGS "-Wl,-rpath,${DEP_INSTALL_DIR}/lib -Wl,-rpath,${GCC_LIB_DIR} -Wl,-rpath,${GCC_LIB32_DIR}"
)

# CPack settings
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Print configuration summary
message(STATUS "Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "  C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "  Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "  OpenMP support: ${OpenMP_CXX_FOUND}")
if(OpenMP_CXX_FOUND)
    message(STATUS "  OpenMP CXX flags: ${OpenMP_CXX_FLAGS}")
endif()
