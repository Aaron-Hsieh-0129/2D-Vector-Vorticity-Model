cmake_minimum_required(VERSION 3.0.0)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(VVM2D VERSION 0.1.0)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# set root of location to find PETSc's pkg-config
set(PETSC $ENV{PETSC_DIR})
set(ENV{PKG_CONFIG_PATH} ${PETSC}/lib/pkgconfig)

# Remove the lines below if you do not wish to have PETSc determine the compilers
execute_process ( COMMAND pkg-config PETSc --variable=ccompiler COMMAND tr -d '\n' OUTPUT_VARIABLE C_COMPILER)
SET(CMAKE_C_COMPILER ${C_COMPILER})
execute_process ( COMMAND pkg-config PETSc --variable=cxxcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE CXX_COMPILER)
if (CXX_COMPILER)
  SET(CMAKE_CXX_COMPILER ${CXX_COMPILER})
endif (CXX_COMPILER)
execute_process ( COMMAND pkg-config PETSc --variable=fcompiler COMMAND tr -d '\n' OUTPUT_VARIABLE FORTRAN_COMPILER)
if (FORTRAN_COMPILER)
  SET(CMAKE_Fortran_COMPILER ${FORTRAN_COMPILER})
  enable_language(Fortran)
endif (FORTRAN_COMPILER)


find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)

set(CMAKE_EXPORT_COMPILE_COMMANDS True)


include_directories(
    include
    /home/Aaron/local/include
    /home/Aaron/local/petsc/include
)

file(GLOB_RECURSE SRC src/*.cpp src/*.hpp)
add_executable(vvm2d ${SRC})

target_link_libraries(vvm2d PkgConfig::PETSC)

find_library(libncPath netcdf "/home/Aaron/local/lib")
find_library(libncxxPath netcdf_c++4 "/home/Aaron/local/lib")
find_library(libpetscPath petsc "/home/Aaron/local/petsc/lib")
target_link_libraries(vvm2d ${libncPath} ${libncxxPath} ${petsc})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
