Tutorials for Installing C/C++ compiler and libraries
======================================================


- Create these directories to put your installed libraries. You can install them anywhere you want. Just keep in mind that you need to add the library paths into ``~/.bashrc`` or ``~/.zshrc`` in this form.

.. code-block:: shell

    mkdir $HOME/local
    mkdir $HOME/local/bin
    mkdir $HOME/local/include
    mkdir $HOME/local/lib

    ################ In ~/.zshrc or ~/.bashrc ##############
    export PATH="$PATH:$HOME/local/bin"
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/local/lib"
    export LIBRARY_PATH="$LIBRARY_PATH:$HOME/local/lib"
    export CPATH="$CPATH:$HOME/local/include"
    #########################################################

Install compiler
-------------------

- If you already have compilers, you can skip this procedure. But the reason I did it here is because I want to have compilers that I have the permission to change without the root privilege.
- Note that if you want to use CUDA, the compiler version for gcc should be less than gcc-10.
.. code-block:: shell

    wget https://ftp.mirrorservice.org/sites/sourceware.org/pub/gcc/releases/gcc-13.2.0/gcc-13.2.0.tar.gz
    tar xf gcc-13.2.0.tar.gz
    mkdir gcc13
    ################ In ~/.zshrc or ~/.bashrc ##############
    unset LIBRARY_PATH LD_LIBRARY_PATH
    #########################################################
    cd gcc-13.2.0

    ./contrib/download_prerequisites

    ./configure --prefix=$HOME/gcc13 --enable-languages=c,c++,fortran --disable-multilib

    make -j64
    make check
    make install

In ``~/.zshrc``

.. code-block:: shell

    ################ In ~/.zshrc or ~/.bashrc ##############
    export PATH=${HOME}/gcc13/bin:${HOME}/gcc13/lib64:$PATH
    export LD_LIBRARY_PATH=${HOME}/gcc13/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=${HOME}/gcc13/lib64:$LD_LIBRARY_PATH
    export LIBRARY_PATH=${HOME}/gcc13/lib:$LIBRARY_PATH
    export LIBRARY_PATH=${HOME}/gcc13/lib64:$LIBRARY_PATH

    #########################################################

Install mpich (Not necessary)
------------------------------
- CUDA can also be installed and link to mpich. You should specify it in the ./configure command.
  
.. code-block:: shell

    wget https://www.mpich.org/static/downloads/4.2.0/mpich-4.2.0.tar.gz
    tar xf mpich-4.2.0.tar.gz
    mkdir mpich4
    cd mpich-4.2.0
    ./configure -prefix=/home/Aaron/mpich4
    make -j64
    make check
    make install
    ################ In ~/.zshrc or ~/.bashrc ##############
    export PATH=/home/Aaron/mpich4/bin:$PATH
    #########################################################


Install cmake 
----------------

.. code-block:: shell

    git clone https://github.com/Kitware/CMake.git
    cd CMake
    git checkout tags/v3.29.1
    ./bootstrap --prefix=$HOME/local

If you encounter some errors that need to resolve by installing some system-level libraries such as openGL or openSSH, you can install cmake from miniconda, pip or install the missing libraries by `sudo apt-get install` from permission of the system administrator.
Building from source is possible and please refer to the official website for more information.


The procedure to install netcdf-cxx
-------------------------------------

Install hdf5
~~~~~~~~~~~~~

.. code-block:: shell

    git clone https://github.com/HDFGroup/hdf5.git
    cd hdf5
    git checkout tags/hdf5-1_8_14
    # ./configure --prefix=$HOME/local --enable-fortran --enable-cxx --enable-parallel --enable-unsupported (error for parallel)
    
    # Use this to enable Fortran
    # ./configure --prefix=$HOME/local --enable-fortran --enable-cxx

    # Use this to enable only C++
    ./configure --prefix=$HOME/local --enable-cxx
    make -j32
    make check
    make install

Install netcdf-c
~~~~~~~~~~~~~~~~~

.. code-block:: shell

    git clone https://github.com/Unidata/netcdf-c.git
    cd netcdf-c
    git checkout tags/v4.3.3.1
    ./configure --prefix=$HOME/local --enable-netcdf-4
    make
    make check
    make install

Noted that if you encounter the error like this `configure: error: Cannot find m4 utility. Install m4 and try again.`
You can install m4 from this website https://ftp.gnu.org/gnu/m4/

.. code-block:: shell

    wget https://ftp.gnu.org/gnu/m4/m4-1.4.1.tar.gz
    tar xf m4-1.4.1.tar.gz
    cd m4-1.4.1
    ./configure --prefix=$HOME/local
    make -j16
    make check
    make install


Install netcdf-cxx
~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    git clone https://github.com/Unidata/netcdf-cxx4.git
    cd netcdf-cxx4
    git checkout v4.2.1
    ./configure --prefix=$HOME/local
    make -j16
    make check
    make install

Install petsc (Not necessary)
----------------

.. code-block:: shell

    git clone -b release https://gitlab.com/petsc/petsc.git petsc
    cd ~/petsc
    ./configure COPTFLAGS="-g -O3" --prefix=${petsc_prefix} --with-openmp=1 --with-cuda=1 --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpifort --download-f2cblaslapack=1 --download-triangle=1

.. code-block:: shell

    ################ In ~/.zshrc or ~/.bashrc ##############
    export PATH=/home/Aaron/gcc13/bin:/home/Aaron/gcc13/lib64:$PATH
    unset LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=/home/Aaron/gcc13/lib/

    export PATH=$PATH:/install/bin
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:/install/include
    export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/install/include
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/install/lib64
    export LIBRARY_PATH=$LIBRARY_PATH:/install/lib

    export PATH="$PATH:$HOME/local/bin"
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/local/lib"
    export LIBRARY_PATH="$LIBRARY_PATH:$HOME/local/lib"
    export CPATH="$CPATH:$HOME/local/include"

    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/gcc13/lib64"
    export LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/gcc13/lib64"

    export PATH=/home/Aaron/mpich4/bin:$PATH

    export petsc_prefix=$HOME/local/petsc
    export PETSC_DIR=/home/Aaron/petsc
    PETSC_ARCH="linux-opt"
    export PYTHONPATH=${petsc_prefix}/lib

    export C_INCLUDE_PATH=$C_INCLUDE_PATH:${HOME}/local/petsc/include
    export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${HOME}/local/include

    export C_INCLUDE_PATH=$C_INCLUDE_PATH:${HOME}/mpich4/include
    export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:${HOME}/mpich4/include
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/mpich4/lib
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${HOME}/mpich4/lib64
    #########################################################

