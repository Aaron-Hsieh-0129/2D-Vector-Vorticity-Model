# Two-dimension Vector Vorticity Model
This is a 2D cloud-resolving model based on vorticity equation

## How to Use
1. Clone the project using
    ```
    git clone https://github.com/Aaron-Hsieh-0129/2D-Vector-Vorticity-Model.git
    ```

2. Using miniconda (or anaconda) to install netcdf_cxx4 library (using the command under terminal)
    ```
    conda install -c conda-forge libnetcdf
    conda install -c conda-forge netcdf-cxx4=4.3.1
    ```

3. Find the installed netcdf library under miniconda and use the paths to modify the paths into CMakeLists.txt. Usually, it will be at ``` <your_path_to_miniconda>/lib/netcdf_c++4 ```

    ```
    find_library(libncxxPath netcdf_c++4 "<your path>")
    ```

4. You are able to run the model by running the command under the project folder
    ```
    sh run.sh
    ```

5. You can change the model settings by  changing the macro flags in the ``` ./src/define.hpp ```

### If you cannot solve the netcdf_cxx4 problem
1. Turn off the flag ```OUTPUTNC``` and turn on the flag ```OUTPUTTXT``` in ``` ./src/define.hpp ```
2. Create a folder called bin under project root by ```mkdir bin```
3. Create a file called Makefile under the project and the content will be the following
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
4. Using ```make``` under project root to compile the project and you will the the execution file at ```./bin/vvm2d``` .
5. You will be able to use the model's output by the txt files


## Develope Status
- This is an ongoing work. The linerized governing equations set is tested and the non-linearized version is not correct now. If you want to use this model, you need to make sure you the macro in ```src/Define.hpp``` is set to be linearized. (2024.02.18)