cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(VVM2D VERSION 0.1.0 LANGUAGES Fortran CXX)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set the C/C++ compiler to gcc/g++ (can be overridden by PETSc if enabled)
set(CMAKE_C_COMPILER "/home/Aaron/gcc13/bin/gcc")
set(CMAKE_CXX_COMPILER "/home/Aaron/gcc13/bin/g++")
set(CMAKE_Fortran_COMPILER "/home/Aaron/gcc13/bin/gfortran")

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fPIC -fdefault-real-8")

# Option to enable or disable PETSc
option(USE_PETSC "Enable PETSc support" OFF)

if(USE_PETSC)
    # Set root location to find PETSc's pkg-config
    set(PETSC $ENV{PETSC_DIR})
    set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

    # PETSc determines compilers if available
    execute_process(COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
    if(C_COMPILER)
        set(CMAKE_C_COMPILER ${C_COMPILER})
    endif()
    
    execute_process(COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
    if(CXX_COMPILER)
        set(CMAKE_CXX_COMPILER ${CXX_COMPILER})
    endif()

    execute_process(COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
    if(FORTRAN_COMPILER)
        set(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
    endif()

    # Find PETSc using pkg-config
    find_package(PkgConfig REQUIRED)
    pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
endif()

# Enable exporting compile commands
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

# Find and enable OpenMP if available
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Include directories
include_directories(
    include
    /home/Aaron/local/include
)

# Add PETSc include directories if PETSc is enabled
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

# Link libraries to vvm2d
target_link_libraries(vvm2d PRIVATE 
    p3  # Link the P3 Fortran library
    gfortran  # Fortran runtime
    ${libncPath}
    ${libncxxPath}
)

# Link PETSc if enabled
if(USE_PETSC)
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
