# Set minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Ensure NVHPC compilers are used for C/C++
set(CMAKE_C_COMPILER nvc CACHE STRING "Set C compiler")
set(CMAKE_CXX_COMPILER nvc++ CACHE STRING "Set C++ compiler")
# Use gfortran for Fortran
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "Set Fortran compiler")

# Disable dependency file generation (-MF issue) for NVHPC
if(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    set(CMAKE_CXX_DEPENDS "")
    set(CMAKE_C_DEPENDS "")
    set(CMAKE_DEPFILE_FLAGS_CXX "")
    set(CMAKE_DEPFILE_FLAGS_C "")
endif()

# Define project after compiler settings
project(VVM2D VERSION 0.1.0 LANGUAGES CXX Fortran)

# Allow user to specify NVHPC and GCC installation paths
set(NVHPC_HOME "" CACHE STRING "Path to NVHPC installation")
set(GCC_HOME "" CACHE STRING "Path to GCC installation")

# Set default paths if not provided
if(NVHPC_HOME STREQUAL "")
    set(NVHPC_HOME "/home/Aaron/nvhpc/Linux_x86_64/21.3" CACHE STRING "Path to NVHPC installation" FORCE)
endif()
if(GCC_HOME STREQUAL "")
    set(GCC_HOME "/home/Aaron/gcc9" CACHE STRING "Path to GCC installation" FORCE)
endif()

set(NVHPC_LIB_DIR "${NVHPC_HOME}/compilers/lib")
set(CUDA_INCLUDE_DIR "${NVHPC_HOME}/cuda/include")
set(CUDA_LIB_DIR "${NVHPC_HOME}/cuda/lib64")
set(GCC9_LIB_DIR "${GCC_HOME}/lib64")
set(GCC9_LIB32_DIR "${GCC_HOME}/lib")


# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set HOME_DIR
set(HOME_DIR $ENV{HOME})

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

# Find MPI (force NVHPC OpenMPI for C/C++)
find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${NVHPC_HOME}/comm_libs/mpi/bin/mpicc)
set(CMAKE_CXX_COMPILER ${NVHPC_HOME}/comm_libs/mpi/bin/mpicxx)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 ${MPI_CXX_COMPILE_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -fdefault-real-8")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

# Add OpenACC flags for GPU support
option(USE_OPENACC "Enable OpenACC GPU support" ON)
if(USE_OPENACC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc -gpu=ccall -Minfo=accel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -gpu=ccall -Minfo=accel")
endif()

# Set RPATH usage
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Include directories
include_directories(
    include
    ${DEP_INSTALL_DIR}/include
    ${MPI_INCLUDE_PATH}
)

# HDF5
find_library(HDF5_LIB NAMES hdf5 hdf5_serial HINTS ${DEP_INSTALL_DIR}/lib NO_DEFAULT_PATH)
if(NOT HDF5_LIB)
    message(STATUS "HDF5 not found in ${DEP_INSTALL_DIR}, installing...")
    ExternalProject_Add(
        hdf5_ext
        URL            https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.gz
        PREFIX         ${CMAKE_SOURCE_DIR}/_deps/hdf5
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
            CC=${CMAKE_C_COMPILER}
            CXX=${CMAKE_CXX_COMPILER}
            <SOURCE_DIR>/configure
            --prefix=${DEP_INSTALL_DIR}
            --enable-cxx
            --disable-shared
            --enable-static
        BUILD_COMMAND  make -j${N_CORES}
        INSTALL_COMMAND make install
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        BUILD_ALWAYS   FALSE
    )
    set(HDF5_NEEDS_INSTALL TRUE)
    set(HDF5_LIB ${DEP_INSTALL_DIR}/lib)
else()
    message(STATUS "HDF5 found in ${DEP_INSTALL_DIR}/lib")
endif()

# NetCDF-C (depends on HDF5)
find_library(NETCDF_LIB NAMES netcdf HINTS ${DEP_INSTALL_DIR}/lib NO_DEFAULT_PATH)
if(NOT NETCDF_LIB)
    message(STATUS "NetCDF-C not found in ${DEP_INSTALL_DIR}, installing...")
    ExternalProject_Add(
        netcdf-c_ext
        GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
        GIT_TAG        v4.9.2
        GIT_SHALLOW    TRUE
        PREFIX         ${CMAKE_SOURCE_DIR}/_deps/netcdf-c
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${DEP_INSTALL_DIR} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} FC=${CMAKE_Fortran_COMPILER} --disable-dap --disable-byterange --with-hdf5=${DEP_INSTALL_DIR}
        BUILD_COMMAND  make -j${N_CORES}
        INSTALL_COMMAND make install
        DEPENDS        hdf5_ext
        UPDATE_COMMAND ""
        BUILD_ALWAYS   FALSE
    )
    set(NETCDF_C_NEEDS_INSTALL TRUE)
    set(NETCDF_LIB ${DEP_INSTALL_DIR}/lib)
else()
    message(STATUS "NetCDF-C found in ${DEP_INSTALL_DIR}/lib")
endif()

# NetCDF-C++4 (depends on NetCDF-C and HDF5)
find_library(NETCDF_CXX4_LIB NAMES netcdf_c++4 HINTS ${DEP_INSTALL_DIR}/lib NO_DEFAULT_PATH)
if(NOT NETCDF_CXX4_LIB)
    message(STATUS "NetCDF-C++4 not found in ${DEP_INSTALL_DIR}, installing...")
    ExternalProject_Add(
        netcdf-cxx4_ext
        GIT_REPOSITORY https://github.com/Unidata/netcdf-cxx4.git
        GIT_TAG        v4.3.1
        GIT_SHALLOW    TRUE
        PREFIX         ${CMAKE_SOURCE_DIR}/_deps/netcdf-cxx4
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${DEP_INSTALL_DIR} --with-netcdf=${DEP_INSTALL_DIR} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CPPFLAGS=-I${DEP_INSTALL_DIR}/include LDFLAGS=-L${DEP_INSTALL_DIR}/lib
        BUILD_COMMAND  make -j${N_CORES}
        INSTALL_COMMAND make install
        DEPENDS        netcdf-c_ext hdf5_ext
        UPDATE_COMMAND ""
        BUILD_ALWAYS   FALSE
    )
    set(NETCDF_CXX4_NEEDS_INSTALL TRUE)
    set(NETCDF_CXX4_LIB ${DEP_INSTALL_DIR}/lib)
else()
    message(STATUS "NetCDF-C++4 found in ${DEP_INSTALL_DIR}/lib")
endif()

# AMGX (version 2.4.0)
find_library(AMGXSH_LIBRARY NAMES amgxsh HINTS ${CMAKE_SOURCE_DIR}/_deps/amgx/src/lib NO_DEFAULT_PATH)
if(NOT AMGXSH_LIBRARY)
    message(STATUS "AMGX not found in ${DEP_INSTALL_DIR}, installing...")
    set(AMGX_NEEDS_INSTALL TRUE)
    ExternalProject_Add(
        amgx_ext
        GIT_REPOSITORY https://github.com/NVIDIA/AMGX.git
        GIT_TAG        v2.4.0
        GIT_SHALLOW    TRUE
        PREFIX         ${CMAKE_SOURCE_DIR}/_deps/amgx
        PATCH_COMMAND  find <SOURCE_DIR> -type f -name "CMakeLists.txt" -exec sed -i "s/-Wno-terminate//g" {} + &&
                   find <SOURCE_DIR> -type f -name "*.cmake" -o -name "CMakeLists.txt" -exec sed -i "s/-mp/-fopenmp/g" {} +
        CMAKE_ARGS     -DCMAKE_CXX_COMPILER=${GCC_HOME}/bin/g++
                       -DCMAKE_C_COMPILER=${GCC_HOME}/bin/gcc
                       -DCMAKE_CUDA_COMPILER=${NVHPC_HOME}/compilers/bin/nvcc
                       -DCMAKE_CUDA_HOST_COMPILER=${GCC_HOME}/bin/g++
                       -DCUDA_TOOLKIT_ROOT_DIR=${NVHPC_HOME}/cuda/11.2
                       -DCUDA_INCLUDE_DIRS=${NVHPC_HOME}/cuda/11.2/include
                       -DCUDA_CUDART_LIBRARY=${NVHPC_HOME}/cuda/11.2/lib64/libcudart.so
                       -DCUDA_cublas_LIBRARY=${NVHPC_HOME}/math_libs/11.2/lib64/libcublas.so
                       -DCUDA_cusparse_LIBRARY=${NVHPC_HOME}/math_libs/11.2/lib64/libcusparse.so
                       -DCUDA_nvToolsExt_LIBRARY=${NVHPC_HOME}/cuda/11.2/lib64/libnvToolsExt.so
                       -DMPI_C_COMPILER=${NVHPC_HOME}/comm_libs/openmpi/openmpi-3.1.5/bin/mpicc
                       -DMPI_CXX_COMPILER=${NVHPC_HOME}/comm_libs/openmpi/openmpi-3.1.5/bin/mpicxx
                       -DMPIEXEC=${NVHPC_HOME}/comm_libs/openmpi/openmpi-3.1.5/bin/mpirun
                       -DMPI_C_INCLUDE_PATH=${NVHPC_HOME}/comm_libs/openmpi/openmpi-3.1.5/include
                       -DMPI_CXX_INCLUDE_PATH=${NVHPC_HOME}/comm_libs/openmpi/openmpi-3.1.5/include
                       -DMPI_C_LIBRARIES=${NVHPC_HOME}/comm_libs/openmpi/openmpi-3.1.5/lib/libmpi.so
                       -DMPI_CXX_LIBRARIES=${NVHPC_HOME}/comm_libs/openmpi/openmpi-3.1.5/lib/libmpi.so
                       -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/_deps/install
                       -DCMAKE_CUDA_ARCHITECTURES=60  # Match your GPU (e.g., Pascal)
                       -DCMAKE_CUDA_FLAGS="-arch=sm_60"  # Enforce compute_60
                       -DCMAKE_CXX_STANDARD=14  # Explicitly enforce C++14
                       -DCMAKE_EXE_LINKER_FLAGS="-L${GCC9_LIB_DIR} -Wl,-rpath,${GCC9_LIB_DIR}"
        BUILD_COMMAND  make -j${N_CORES}
        INSTALL_COMMAND make install 
        UPDATE_COMMAND ""
        BUILD_ALWAYS   FALSE
    )

    set(AMGXSH_LIBRARY "${DEP_INSTALL_DIR}/lib/libamgxsh.so")
else()
    message(STATUS "AMGX found in ${CMAKE_SOURCE_DIR}/amgx/src/lib")
endif()

# Option to enable or disable PETSc
option(USE_PETSC "Enable PETSc support" OFF)
if(USE_PETSC)
    message(STATUS "PETSc support enabled")
    find_library(PETSC_LIB NAMES petsc HINTS ${DEP_INSTALL_DIR}/lib NO_DEFAULT_PATH)
    if(NOT PETSC_LIB)
        message(STATUS "PETSc not found in ${DEP_INSTALL_DIR}, installing...")
        ExternalProject_Add(
            petsc_ext
            GIT_REPOSITORY https://gitlab.com/petsc/petsc.git
            GIT_TAG        v3.20.5
            GIT_SHALLOW    TRUE
            PREFIX         ${CMAKE_SOURCE_DIR}/_deps/petsc
            CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${DEP_INSTALL_DIR} --with-cc=${CMAKE_C_COMPILER} --with-cxx=${CMAKE_CXX_COMPILER} --with-fc=${CMAKE_Fortran_COMPILER} --with-mpi=1
            BUILD_COMMAND  make -j${N_CORES}
            INSTALL_COMMAND make install
            UPDATE_COMMAND ""
            BUILD_ALWAYS   FALSE
        )
        set(PETSC_NEEDS_INSTALL TRUE)
        set(PETSC_LIB ${DEP_INSTALL_DIR}/lib/libpetsc.so)
    else()
        message(STATUS "PETSc found in ${DEP_INSTALL_DIR}/lib")
    endif()
endif()

# CUDA settings

include_directories(
    include
    ${MPI_INCLUDE_PATH}
    ${CUDA_INCLUDE_DIR}
    ${DEP_INSTALL_DIR}/include
)

# Create the P3 Fortran shared library
add_library(p3 SHARED src/microphy_p3.f90)
set_target_properties(p3 PROPERTIES
    OUTPUT_NAME "p3"
    Fortran_PREPROCESS ON
)

# Source files for VVM2D
file(GLOB_RECURSE SRC src/*.cpp src/*.hpp)
add_executable(vvm2d ${SRC})

# Link Fortran library and runtime to the executable
target_link_libraries(vvm2d PRIVATE p3 gfortran)

# Link MPI libraries (C/C++ only)
target_link_libraries(vvm2d PRIVATE ${MPI_CXX_LIBRARIES})

# Link libraries
find_library(libncPath netcdf "${NETCDF_LIB}")
find_library(libncxxPath netcdf_c++4 "${NETCDF_CXX4_LIB}")
find_library(cudartPath cudart ${CUDA_LIB_DIR})

target_link_libraries(vvm2d PRIVATE
    ${libncPath}
    ${libncxxPath}
    ${AMGXSH_LIBRARY}
    ${cudartPath}
    stdc++
)

# Link PETSc if enabled
if(USE_PETSC)
    find_library(libpetscPath petsc "${PETSC_LIB}")
    target_link_libraries(vvm2d PRIVATE libpetscPath)
endif()

# Add dependencies if needed
if(HDF5_NEEDS_INSTALL)
    add_dependencies(vvm2d hdf5_ext)
endif()
if(NETCDF_C_NEEDS_INSTALL)
    add_dependencies(vvm2d netcdf-c_ext)
endif()
if(NETCDF_CXX4_NEEDS_INSTALL)
    add_dependencies(vvm2d netcdf-cxx4_ext)
endif()
if(AMGX_NEEDS_INSTALL)
    add_dependencies(vvm2d amgx_ext)
endif()
if(USE_PETSC AND PETSC_NEEDS_INSTALL)
    add_dependencies(vvm2d petsc_ext)
endif()

# Set runtime library paths (-rpath)
set_target_properties(vvm2d PROPERTIES
    LINK_FLAGS "-Wl,-rpath,${GCC9_LIB_DIR} -Wl,-rpath,${GCC9_LIB32_DIR} -Wl,-rpath,${CUDA_LIB_DIR} -Wl,-rpath,${DEP_INSTALL_DIR}/lib"
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
message(STATUS "  PETSc support: ${USE_PETSC}")
message(STATUS "  OpenACC support: ${USE_OPENACC}")
