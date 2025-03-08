cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(VVM2D VERSION 0.1.0 LANGUAGES Fortran CXX)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set the C/C++ compiler to nvc++ (NVIDIA HPC SDK)
set(CMAKE_C_COMPILER "nvc")
set(CMAKE_CXX_COMPILER "nvc++")

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Add OpenACC flags for GPU support
option(USE_OPENACC "Enable OpenACC GPU support" ON)
if(USE_OPENACC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc -gpu=ccall -Minfo=accel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -gpu=ccall -Minfo=accel")
endif()

# Disable automatic dependency tracking for NVIDIA compilers
set(CMAKE_C_DEPENDS_USE_COMPILER FALSE)
set(CMAKE_CXX_DEPENDS_USE_COMPILER FALSE)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -fdefault-real-8")

# Option to enable or disable PETSc
option(USE_PETSC "Enable PETSc support" OFF)

# AMGX and CUDA settings
set(AMGX_INCLUDE_DIR "/home/Aaron/amgx_install/include")
set(AMGX_LIB_DIR "/home/Aaron/amgx_install/lib")
set(CUDA_INCLUDE_DIR "/home/Aaron/nvhpc/Linux_x86_64/21.1/cuda/include")
set(CUDA_LIB_DIR "/home/Aaron/nvhpc/Linux_x86_64/21.1/cuda/lib64")
set(GCC9_LIB_DIR "/home/Aaron/gcc9/lib64")
set(GCC9_LIB32_DIR "/home/Aaron/gcc9/lib")  # Added for completeness

# Include directories
include_directories(
    include
    /home/Aaron/local/include
    ${AMGX_INCLUDE_DIR}
    ${CUDA_INCLUDE_DIR}
)

# Add PETSc include directories if enabled
if(USE_PETSC)
    include_directories(/home/Aaron/local/petsc/include)
endif()

# Create the P3 Fortran shared library
add_library(p3 SHARED src/microphy_p3.f90)
set_target_properties(p3 PROPERTIES OUTPUT_NAME "p3")

# Source files for VVM2D
file(GLOB_RECURSE SRC src/*.cpp src/*.hpp)
add_executable(vvm2d ${SRC})

# Find libraries
find_library(libncPath netcdf "/home/Aaron/local/lib")
find_library(libncxxPath netcdf_c++4 "/home/Aaron/local/lib")
find_library(libpetscPath petsc "/home/Aaron/local/petsc/lib")
find_library(AMGXSH_LIBRARY amgxsh ${AMGX_LIB_DIR})
find_library(CUDART_LIBRARY cudart ${CUDA_LIB_DIR})
find_library(GFORTRAN_LIBRARY gfortran "/usr/lib64" HINTS "${GCC9_LIB_DIR}")

# Ensure nvc++ uses GCC 9's libstdc++
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L${GCC9_LIB_DIR} -L${GCC9_LIB32_DIR}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${GCC9_LIB_DIR} -L${GCC9_LIB32_DIR}")

# Link libraries to vvm2d
target_link_libraries(vvm2d PRIVATE 
    p3  # Link the P3 Fortran library
    ${GFORTRAN_LIBRARY}  # Fortran runtime
    ${libncPath}
    ${libncxxPath}
    ${AMGXSH_LIBRARY}  # AMGX shared library
    ${CUDART_LIBRARY}  # CUDA runtime library
    stdc++  # Explicitly link libstdc++
)

# Add runtime library paths (-rpath)
set_target_properties(vvm2d PROPERTIES
    LINK_FLAGS "-Wl,-rpath,${GCC9_LIB_DIR} -Wl,-rpath,${GCC9_LIB32_DIR} -Wl,-rpath,${CUDA_LIB_DIR} -Wl,-rpath,${AMGX_LIB_DIR}"
)

# Link PETSc if enabled
if(USE_PETSC)
    set(ENV{PKG_CONFIG_PATH} $ENV{PETSC_DIR}/lib/pkgconfig)
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
    target_link_libraries(vvm2d PRIVATE 
        PkgConfig::PETSC
        ${libpetscPath}
    )
endif()

# Ensure the P3 library is built before vvm2d
add_dependencies(vvm2d p3)

# CPack settings
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
