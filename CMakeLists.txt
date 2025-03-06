cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(VVM2D VERSION 0.1.0 LANGUAGES Fortran CXX)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Set the C/C++ compiler to pgcc/pgc++
set(CMAKE_C_COMPILER "pgcc")
set(CMAKE_CXX_COMPILER "pgc++")

# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# Add OpenACC flags for GPU support
option(USE_OPENACC "Enable OpenACC GPU support" ON)  # Enabled by default, can be turned off with -DUSE_OPENACC=OFF
if(USE_OPENACC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc -gpu=ccall -Minfo=accel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -gpu=ccall -Minfo=accel")
endif()

# Disable automatic dependency tracking for PGI compilers
set(CMAKE_C_DEPENDS_USE_COMPILER FALSE)
set(CMAKE_CXX_DEPENDS_USE_COMPILER FALSE)
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

# Enable exporting compile commands (optional, but may not work fully with PGI)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

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
