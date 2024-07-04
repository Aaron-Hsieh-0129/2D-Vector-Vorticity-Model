# Two-dimension Vector Vorticity Model

This is a 2D cloud-resolving model based on vorticity equation

- [Two-dimension Vector Vorticity Model](#two-dimension-vector-vorticity-model)
  - [Prerequisite](#prerequisite)
  - [How to Use](#how-to-use)

## Prerequisite

- C++ compiler (higher than C++11)
- netcdf-cxx4 (hdf5, netcdf-c are needed for netcdf-cxx) [optional]
- PETSc [optional]
- Eigen (this has already be installed in include folder) [optional]

The tutorial for installing netcdf-cxx4 and PETSc can be found [here](./Install_compilers_libraries.md)

- This model will use txt output and Eigen solver for solving Poisson equation in default. However
  - You can turn on the `OUTPUTNC` and turn off `OUTPUTTXT` in `./src/Config.hpp` to use netcdf output. Note that the netcdf doesn't support openMP output so if you want to use openMP, don't turn the `OUTPUTNC` flag on.
  - Turn on `PETSC` in `./src/Config.hpp`, the model will change the Poisson solver package to `PETSc`
  - 
## How to Use

1. Clone the project using

   ```
   git clone https://github.com/Aaron-Hsieh-0129/2D-Vector-Vorticity-Model.git
   ```

2. Install netcdf-cxx, petsc [optional]

   It's a little bit complicated to install libraries for C/C++.

   I will provide a tutorial for installing C/C++ compiler and the libraries in another file [here](./Install_compilers_libraries.md).
   Here, you don't need to have sudo privilege to install anything.

3. Link the installed libraries [optional]

- If you don't need to use PETSc and netcdf, you might need to turn off the link command in CMakeLists.txt.

- You need to change the libraries path (netcdf, petsc) to your own path.
- Change include path in CMakeLists.txt
  ```CMake
  include_directories(
    include
    </path/to/your/petsc>/include
  )
  ```
- Change library link path
  ```
  find_library(libncxxPath netcdf_c++4 "<path to your netcdf_c++4>/lib")
  find_library(libpetscPath petsc "<path to your petsc>/lib")
  ```

1. You are able to run the model by running the command under the project folder

   ```
   sh run.sh
   ```

   or you can use your own command by referencing the command in `run.sh`

2. You can change the model settings by changing the macro flags in the `./src/Config.hpp`
