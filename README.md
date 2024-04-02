# Two-dimension Vector Vorticity Model
This is a 2D cloud-resolving model based on vorticity equation

- [Two-dimension Vector Vorticity Model](#two-dimension-vector-vorticity-model)
  - [Prerequisite](#prerequisite)
  - [How to Use](#how-to-use)
    - [If you cannot solve the netcdf\_cxx4 and petsc installation problem.](#if-you-cannot-solve-the-netcdf_cxx4-and-petsc-installation-problem)

## Prerequisite
- C++ compiler (higher than C++11)
- netcdf-cxx4 (hdf5, netcdf-c are needed for netcdf-cxx)
- PETSc
- Eigen (this has already be installed in include folder)

The tutorial for installing them can be found [here](./Install_compilers_libraries.md)

- This model will use netcdf-cxx and petsc in default. However
  - You can turn off the `OUTPUTNC`  and turn on  `OUTPUTTXT` in `./src/Config.hpp` to use txt output without installing `netcdf-cxx4`.
  - Turn off `PETSC` in `./src/Config.hpp`, the model will change the Poisson solver package to `Eigen`, which means you don't need to install PETSc.

## How to Use
1. Clone the project using
    ```
    git clone https://github.com/Aaron-Hsieh-0129/2D-Vector-Vorticity-Model.git
    ```

2. Install netcdf-cxx, petsc
   
   It's a little bit complicated to install libraries for C/C++.
   
   I will provide a tutorial for installing C/C++ compiler and the libraries in another file [here](./Install_compilers_libraries). 
   Here, you don't need to have sudo privilege to install anything.


3. Link the installed libraries 

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

2. You can change the model settings by  changing the macro flags in the ``` ./src/Config.hpp ```

### If you cannot solve the netcdf_cxx4 and petsc installation problem.
1. Turn off the flag ```OUTPUTNC``` and turn on the flag ```OUTPUTTXT``` in ``` ./src/Config.hpp ```
2. Turn off the flag ```PETSC``` in ``` ./src/Config.hpp ``` and this will in turn make the model use the Eigen C++ Solver
3. Create a folder called bin under project root by ```mkdir bin```
4. Create a file called Makefile under the project and the content will be the following
    ```Makefile
    # CC and CFLAGS are varilables
    CC = g++
    CFLAGS = -Wall -c -std=c++11
    # -c option ask g++ to compile the source files, but do not link.
    OPTFLAGS = -O3

    all	: bin/vvm2d
        @echo -n ""

    bin/vvm2d	: main.o Declare.o Init.o Iteration.o Outputfile.o 
                $(CC) $(OPTFLAGS) main.o Declare.o Init.o Outputfile.o Iteration.o -o bin/vvm2d
    main.o 	   	: src/main.cpp
                $(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
    Declare.o	: src/Declare.cpp
                $(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
    Init.o		: src/Init.cpp
                $(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
    Outputfile.o: src/Outputfile.cpp
                $(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
    Iteration.o	: src/Iteration.cpp
                $(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@
    ```
5. Using ```make``` under project root to compile the project and you will the the execution file at ```./bin/vvm2d``` .
6. You will be able to use the model's output by the txt files
