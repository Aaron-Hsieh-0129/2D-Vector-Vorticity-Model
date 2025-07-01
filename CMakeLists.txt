# Set minimum CMake version
cmake_minimum_required(VERSION 3.10)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()


# Allow user to specify NVHPC and GCC installation paths
set(NVHPC_HOME "/work/aaron900129/nvhpc/Linux_x86_64/21.3" CACHE STRING "Path to NVHPC installation")
set(GCC_HOME "/work/aaron900129/gcc9" CACHE STRING "Path to GCC installation")

# Set compilers explicitly before project()
set(CMAKE_C_COMPILER "${NVHPC_HOME}/comm_libs/mpi/bin/mpicc" CACHE PATH "C compiler" FORCE)
set(CMAKE_CXX_COMPILER "${NVHPC_HOME}/comm_libs/mpi/bin/mpicxx" CACHE PATH "C++ compiler" FORCE)
set(CMAKE_Fortran_COMPILER "${GCC_HOME}/bin/gfortran" CACHE PATH "Fortran compiler" FORCE)
set(CMAKE_CUDA_COMPILER "${NVHPC_HOME}/compilers/bin/nvcc" CACHE PATH "CUDA compiler" FORCE)

# Disable dependency file generation (-MF issue) for NVHPC
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    set(CMAKE_CXX_DEPENDS "")
    set(CMAKE_C_DEPENDS "")
    set(CMAKE_DEPFILE_FLAGS_CXX "")
    set(CMAKE_DEPFILE_FLAGS_C "")
endif()

# Define project after compiler settings
project(VVM2D VERSION 0.1.0 LANGUAGES CXX Fortran)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set HOME_DIR
set(HOME_DIR "/work/aaron900129")

# Find MPI (force NVHPC OpenMPI for C/C++)
# set(ENV{OMPI_CC} nvc)
# set(ENV{OMPI_CXX} nvc++)
# set(ENV{OMPI_FC} "/work/aaron900129/gcc9/bin/gfortran") # Set Fortran compiler for MPI, though not used in Fortran
find_package(MPI REQUIRED)

# Set MPI compilers only for C/C++ (Fortran compiler remains gfortran)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 ${MPI_CXX_COMPILE_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -fdefault-real-8") # Adjusted to gfortran flags
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

# Add OpenACC flags for GPU support
option(USE_OPENACC "Enable OpenACC GPU support" ON)
if(USE_OPENACC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc -gpu=ccall -Minfo=accel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -gpu=ccall -Minfo=accel")
endif()

# Set RPATH usage
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_definitions(-DRTE_USE_CBOOL)

# Option to enable or disable PETSc
option(USE_PETSC "Enable PETSc support" OFF)

if(USE_PETSC)
    # Set root location to find PETSc's pkg-config
    set(PETSC $ENV{PETSC_DIR})
    set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
endif()

# AMGX and CUDA settings
set(AMGX_INCLUDE_DIR "${HOME_DIR}/amgx_install/include")
set(AMGX_LIB_DIR "${HOME_DIR}/amgx_install/lib")
set(CUDA_INCLUDE_DIR "${HOME_DIR}/nvhpc/Linux_x86_64/21.3/cuda/include")
set(CUDA_LIB_DIR "${HOME_DIR}/nvhpc/Linux_x86_64/21.3/cuda/lib64")
set(GCC9_LIB_DIR "${HOME_DIR}/gcc9/lib64")
set(GCC9_LIB32_DIR "${HOME_DIR}/gcc9/lib")  # Added for completeness

# Include directories
include_directories(
    include
    ${HOME_DIR}/local/include
    ${MPI_INCLUDE_PATH}
    ${AMGX_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include
    ${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include_test
)

# Add PETSc include directories if enabled
if(USE_PETSC)
    include_directories(${HOME_DIR}/local/petsc/include)
endif()

# Create the P3 Fortran shared library
add_library(p3 SHARED external/P3-microphysics/microphy_p3.f90)
set_target_properties(p3 PROPERTIES OUTPUT_NAME "p3")

# Source files for VVM2D
file(GLOB_RECURSE SRC src/*.cpp src/*.hpp)
add_executable(vvm2d ${SRC})

# Link Fortran library and runtime to the executable
target_link_libraries(vvm2d PRIVATE p3 gfortran)
add_dependencies(vvm2d p3)

# Link MPI libraries (C/C++ only)
target_link_libraries(vvm2d PRIVATE ${MPI_CXX_LIBRARIES})

# Link PETSc libraries if enabled
if(USE_PETSC)
    target_link_libraries(vvm2d PRIVATE PkgConfig::PETSC)
endif()

# Add rte-rrtmgp-cpp as a subdirectory
add_subdirectory(external/rte-rrtmgp-cpp)

# Force rte-rrtmgp to use only _deps/install headers
set_target_properties(rte_rrtmgp PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${DEP_INSTALL_DIR}/include;${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include"
)


# Find and link other libraries
find_library(libncPath netcdf "${HOME_DIR}/local/lib")
find_library(libncxxPath netcdf_c++4 "${HOME_DIR}/local/lib")
find_library(libpetscPath petsc "${HOME_DIR}/local/petsc/lib")
find_library(AMGXSH_LIBRARY amgxsh ${AMGX_LIB_DIR})
find_library(CUDART_LIBRARY cudart ${CUDA_LIB_DIR})
set(BOOST_SYSTEM_LIB /work/aaron900129/local/lib/libboost_system.so)
set(BOOST_FILESYSTEM_LIB /work/aaron900129/local/lib/libboost_filesystem.so)

# Set runtime library paths (-rpath)
set_target_properties(vvm2d PROPERTIES
    LINK_FLAGS "-Wl,-rpath,${GCC9_LIB_DIR} -Wl,-rpath,${GCC9_LIB32_DIR} -Wl,-rpath,${CUDA_LIB_DIR} -Wl,-rpath,${AMGX_LIB_DIR}"
)

# Link libraries to vvm2d
target_link_libraries(vvm2d PRIVATE
    rte_rrtmgp
    ${libncPath}
    ${libncxxPath}
    ${AMGXSH_LIBRARY}
    ${CUDART_LIBRARY}
    stdc++
    ${BOOST_SYSTEM_LIB}
    ${BOOST_FILESYSTEM_LIB}
)

# Link PETSc libraries conditionally
if(USE_PETSC)
    target_link_libraries(vvm2d PRIVATE ${libpetscPath})
endif()

# CPack settings
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
