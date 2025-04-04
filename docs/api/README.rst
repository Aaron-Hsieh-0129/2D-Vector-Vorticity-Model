Two-dimension Vector Vorticity Model
====================================

This is a 2D cloud-resolving model based on the vorticity equation.


Prerequisite
------------

- C++ compiler (higher than C++11)
- CMake (higher than 3.18.0)
- hdf5 1.8.14 (can be automatically downloaded through CMake)
- netcdf-c 4.3.3.1 (can be automatically downloaded through CMake)
- netcdf-cxx4 4.2.1 (can be automatically downloaded through CMake)
- Eigen (this has already been installed in the include folder)
- PETSc [optional] (PETSc for solving Poisson equation can be used by turn on the flag in `./src/Config.hpp` but in general, Eigen can do it well)
- nvhpc 21.3 [necessary for GPU usage] (If you want to use GPU to run this model, the nvhpc needs to be installed in nvidia website.)
- AMGX 2.4.0 [necessary for GPU usage] (can be automatically downloaded through CMake. Noted that if it has error during download, you might need to specify the nvhpc path and library paths)
If you don't want to use GPU (no nvhpc), please make the the flag for `GPU_POISSON` in `./src/Config.hpp` is not open!!!

It might takes about 20 minutes to install all the dependencies for the first time running.

How to Use
----------

1. Clone the project using:

   .. code-block:: bash
      git clone https://github.com/Aaron-Hsieh-0129/2D-Vector-Vorticity-Model.git 2DVVM

2. You can change most of the model settings in `./vvm_configs.txt`.

3. If you want to use GPU to solve the model, please remember to install nvhpc. You'll need to specify your nvhpc path and gcc (or other supported compiler) path when using CMake to build.

   .. code-block:: bash
      mkdir build && cd build
      # e.g. (CPU) cmake -DGCC_HOME=/home/Aaron/gcc9 ..
      # e.g. (GPU) cmake -DNVHPC_HOME=/home/Aaron/nvhpc/Linux_x86_64/21.3 -DGCC_HOME=/home/Aaron/gcc9 ..
      cmake -DNVHPC_HOME=/path/to/nvhpc -DGCC_HOME=/path/to/gcc ..
      make
      # (CPU) mpirun -np 1 ./vvm2d
      # (GPU) mpirun -np 1 -mca btl_base_warn_component_unused 0 -np 1 -x CUDA_VISIBLE_DEVICES=0 ./vvm2d

    > Please choose the correct cmake file from folder `scripts/`. `CMakeLists_CPU.txt` for CPU case and `CMakeLists_GPU.txt` for GPU.
    > Move the desired one to project root, e.g. `cp scripts/CMakeLists_GPU.txt CMakeLists.txt`. 
    > The CMake will automatically download the dependencies under `_deps` at the first execution.
    > Please make sure the GPU_POISSON flag is commented out in `src/Config.hpp` for non-GPU case.
    > The shell script `run.sh` can be referenced to finish all the process in step 3. 


Any Problems
-------------
- If you encounter any errors when installing the dependencies from cmake, it's recommended to install gcc9.4. There is tutorial for installing gcc in the documentation of 2DVVM.


Documentation
-------------
https://aaron-hsieh-0129.github.io/2D-Vector-Vorticity-Model/index.html
