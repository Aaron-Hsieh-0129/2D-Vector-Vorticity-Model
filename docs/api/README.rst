Two-dimension Vector Vorticity Model
====================================

This is a 2D cloud-resolving model based on the vorticity equation.


Prerequisite
------------

- C++ compiler (higher than C++11)
- CMake (higher than 3.0.0) (You can create your own Makefile by translating the CMakefile.txt if you don't want to use CMake)
- netcdf-cxx4 (hdf5, netcdf-c are needed for netcdf-cxx) [optional]
- PETSc [optional]
- Eigen (this has already been installed in the include folder) [optional]

The tutorial for installing netcdf-cxx4 and PETSc can be found `here <./api/install_compilers_libraries.html>`_ (If the link doesn't work, find this in the 2DVVM documentation)

- This model will use txt output and Eigen solver for solving Poisson equation by default. However,
  - You can turn on the `OUTPUTNC` and turn off `OUTPUTTXT` in `./src/Config.hpp` to use netcdf output. Note that netcdf doesn't support openMP output so if you want to use openMP, don't turn the `OUTPUTNC` flag on.
  - Turn on `PETSC` in `./src/Config.hpp`, the model will change the Poisson solver package to `PETSc`.

How to Use
----------

1. Clone the project using:

   .. code-block:: bash

      git clone https://github.com/Aaron-Hsieh-0129/2D-Vector-Vorticity-Model.git

2. You can change the model settings by changing the macro flags in the `./src/Config.hpp` and also the configuration object in the `./src/main.cpp`. In general, you only need to modify the configuration in `./src/main.cpp`. But if you want to change the numerical methods or the running cases, you might need to modify the flags in the `./src/Config.hpp`.

3. You are able to run the model by running the command under the project folder:

   .. code-block:: bash

      sh run.sh

   or you can use your own command by referencing the command in `run.sh`.

Optional for NetCDF output and PETSc Solver
-------------------------------------------

1. Install netcdf-cxx, petsc.

   It's a little bit complicated to install libraries for C/C++. I will provide a tutorial for installing C/C++ compiler and the libraries in another file `here <./api/install_compilers_libraries.html>`_ (If the link doesn't work, find this in the 2DVVM documentation). Here, you don't need to have sudo privilege to install anything.

2. Link the installed libraries.

- If you don't need to use PETSc and netcdf, you might need to turn off the link command in CMakeLists.txt.
- You need to change the libraries path (netcdf, petsc) to your own path.
- Change include path in CMakeLists.txt:

  .. code-block:: cmake

    include_directories(
      include
      </path/to/your/petsc>/include
    )

- Change library link path:

  .. code-block:: cmake

    find_library(libncxxPath netcdf_c++4 "<path to your netcdf_c++4>/lib")
    find_library(libpetscPath petsc "<path to your petsc>/lib")


3. Change the flags in `./src/Config.hpp`:

- Turn on `OUTPUTNC` and turn off `OUTPUTTXT` to use netcdf output.
- Turn on `PETSC` to use PETSc solver.


Documentation
-------------
https://aaron-hsieh-0129.github.io/2D-Vector-Vorticity-Model/index.html