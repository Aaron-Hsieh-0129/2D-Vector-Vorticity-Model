cmake_minimum_required(VERSION 3.10)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Define project
project(VVM2D VERSION 0.1.0 LANGUAGES CXX Fortran)

# Allow user to specify GCC installation path
set(GCC_HOME "" CACHE STRING "Path to GCC installation")
if(GCC_HOME STREQUAL "")
    set(GCC_HOME "/home/Aaron/gcc9" CACHE STRING "Path to GCC installation" FORCE)
endif()
set(GCC_LIB_DIR "${GCC_HOME}/lib64")
set(GCC_LIB32_DIR "${GCC_HOME}/lib")

# Set GCC compilers (can be overridden by PETSc if enabled)
set(CMAKE_C_COMPILER "${GCC_HOME}/bin/gcc" CACHE STRING "Set C compiler")
set(CMAKE_CXX_COMPILER "${GCC_HOME}/bin/g++" CACHE STRING "Set C++ compiler")
set(CMAKE_Fortran_COMPILER "${GCC_HOME}/bin/gfortran" CACHE STRING "Set Fortran compiler")

# Find MPI (uncomment if needed)
# find_package(MPI REQUIRED)
# set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
# set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

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

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 ${MPI_CXX_COMPILE_FLAGS} -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -fdefault-real-8 -ffree-line-length-0")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")

# Add OpenMP flags
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
if(OpenMP_Fortran_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# Set RPATH usage
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_definitions(-DRTE_USE_CBOOL)

# Include directories
include_directories(
    include
    /home/Aaron/local/include
    ${DEP_INSTALL_DIR}/include
    ${MPI_INCLUDE_PATH}
    ${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include
    ${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include_test
)

# Add rte-rrtmgp-cpp as a subdirectory
add_subdirectory(external/rte-rrtmgp-cpp)

# Force rte-rrtmgp to use only _deps/install headers
set_target_properties(rte_rrtmgp PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "/home/Aaron/local/include;${CMAKE_SOURCE_DIR}/external/rte-rrtmgp-cpp/include"
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

# Link MPI libraries (uncomment if MPI is enabled)
# target_link_libraries(vvm2d PRIVATE ${MPI_CXX_LIBRARIES})

# Link OpenMP libraries
if(OpenMP_CXX_FOUND)
    target_link_libraries(vvm2d PRIVATE OpenMP::OpenMP_CXX)
endif()

# Link libraries
find_library(libncPath netcdf "/home/Aaron/local/lib")
find_library(libncxxPath netcdf_c++4 "/home/Aaron/local/lib")

target_link_libraries(vvm2d PRIVATE
    rte_rrtmgp
    ${libncPath}
    ${libncxxPath}
    stdc++
)

# Link PETSc if enabled
option(USE_PETSC "Enable PETSc support" OFF)
if(USE_PETSC)
    message(STATUS "PETSc support enabled")
    find_library(libpetscPath petsc "${DEP_INSTALL_DIR}/lib" NO_DEFAULT_PATH)
    if(NOT libpetscPath)
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
        set(libpetscPath ${DEP_INSTALL_DIR}/lib/libpetsc.so)
    else()
        message(STATUS "PETSc found in ${DEP_INSTALL_DIR}/lib")
    endif()
    target_link_libraries(vvm2d PRIVATE ${libpetscPath})
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
if(USE_PETSC AND PETSC_NEEDS_INSTALL)
    add_dependencies(vvm2d petsc_ext)
endif()

# Set runtime library paths (-rpath)
set_target_properties(vvm2d PROPERTIES
    LINK_FLAGS "-Wl,-rpath,${GCC_LIB_DIR} -Wl,-rpath,${GCC_LIB32_DIR} -Wl,-rpath,${DEP_INSTALL_DIR}/lib"
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
message(STATUS "  OpenMP support: ${OpenMP_CXX_FOUND}")
if(OpenMP_CXX_FOUND)
    message(STATUS "  OpenMP CXX flags: ${OpenMP_CXX_FLAGS}")
endif()
