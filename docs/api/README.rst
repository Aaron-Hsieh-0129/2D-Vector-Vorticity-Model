Two-Dimension Vector Vorticity Model
====================================

This is a 2D cloud-resolving model based on the vorticity equation.

Prerequisites
-------------

The following dependencies are required:

- **C++ compiler** (version >= C++11)
- **CMake** (version >= 3.18.0)
- **hdf5** 1.8.14 *(can be automatically downloaded via CMake)*
- **netcdf-c** 4.3.3.1 *(can be automatically downloaded via CMake)*
- **netcdf-cxx4** 4.2.1 *(can be automatically downloaded via CMake)*
- **Eigen** *(pre-installed in the include folder)*
- **PETSc** *(optional)* - Required for solving the Poisson equation if enabled in `./src/Config.hpp`. However, Eigen generally suffices.
- **nvhpc** 21.3 *(required for GPU usage)* - Installable from the NVIDIA website.
- **AMGX** 2.4.0 *(required for GPU usage, can be automatically downloaded via CMake)*

  .. note::
     If you encounter errors while downloading AMGX, you may need to specify the NVHPC path and library paths manually.

  .. warning::
     If you do **not** want to use GPU, ensure that the `GPU_POISSON` flag in `./src/Config.hpp` is **disabled**.

Installing dependencies for the first time may take approximately **20 minutes**.

How to Use
----------

1. Clone the project repository:

   .. code-block:: bash

      git clone https://github.com/Aaron-Hsieh-0129/2D-Vector-Vorticity-Model.git 2DVVM

2. Modify model settings as needed in `./vvm_configs.txt`.

3. **Building the Project:**

   If you want to use **GPU**, install NVHPC and specify paths for NVHPC and GCC when configuring with CMake.

   .. code-block:: bash

      mkdir build && cd build
      # Example (CPU): Ensure the `GPU_POISSON` flag is **commented out** in `src/Config.hpp`
      cmake -DGCC_HOME=/home/Aaron/gcc9 ..
      # Example (GPU):
      cmake -DNVHPC_HOME=/home/Aaron/nvhpc/Linux_x86_64/21.3 -DGCC_HOME=/home/Aaron/gcc9 ..
      make
      
      # Running the model:
      # CPU Execution:
      ./vvm2d
      # GPU Execution:
      mpirun -np 1 -mca btl_base_warn_component_unused 0 -np 1 -x CUDA_VISIBLE_DEVICES=0 ./vvm2d

   .. tip::
      - Choose the correct CMake file from the `scripts/` directory:
        - `CMakeLists_CPU.txt` for CPU execution.
        - `CMakeLists_GPU.txt` for GPU execution.
      - Copy the appropriate file to the project root:

        .. code-block:: bash

           cp scripts/CMakeLists_GPU.txt CMakeLists.txt

      - On first execution, CMake will download dependencies to the `_deps/` directory automatically.
      - Ensure the `GPU_POISSON` flag is **commented out** in `src/Config.hpp` if using **CPU only**.
      - You can use the `run.sh` script to automate these steps.

Troubleshooting
---------------

- If dependency installation fails via CMake, installing **GCC 9.4** is recommended. 
- Refer to the official documentation for GCC installation steps.

Documentation
-------------

For more details, visit:

`2D Vector Vorticity Model Documentation <https://aaron-hsieh-0129.github.io/2D-Vector-Vorticity-Model/index.html>`_

