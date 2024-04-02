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

.. code-block:: shell

    wget https://ftp.mirrorservice.org/sites/sourceware.org/pub/gcc/releases/gcc-13.2.0/gcc-13.2.0.tar.gz
    tar xf gcc-13.2.0.tar.gz
    mkdir gcc13
    ################ In ~/.zshrc or ~/.bashrc ##############
    unset LIBRARY_PATH LD_LIBRARY_PATH
    #########################################################
    cd gcc-13.2.0

    ./contrib/download_prerequisites

    ./configure --prefix=$HOME/gcc13 --enable-languages=c,c++,fortran,go --disable-multilib

    make -j 60
    make check
    make install

In ``~/.zshrc``

.. code-block:: shell

    ################ In ~/.zshrc or ~/.bashrc ##############
    export PATH=/path/to/software/gcc9/bin:/path/to/software/gcc9/lib64:$PATH
    export LD_LIBRARY_PATH=/path/to/software/gcc9/lib/:$LD_LIBRARY_PATH
    export PATH=$PATH:/install/bin
    export C_INCLUDE_PATH=$C_INCLUDE_PATH:/install/include
    export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/install/include
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/install/lib64
    export LIBRARY_PATH=$LIBRARY_PATH:/install/lib

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/gcc13/lib64
    export LIBRARY_PATH="$LIBRARY_PATH:$HOME/gcc13/lib64"
    #########################################################

Install mpich
----------------

.. code-block:: shell

    wget https://www.mpich.org/static/downloads/4.2.0/mpich-4.2.0.tar.gz
    tar xf mpich-4.2.0.tar.gz
    mkdir mpich4
    cd mpich-4.2.0
    ./configure -prefix=/home/Aaron/mpich4
    make -j 64
    make check
    make install
    ################ In ~/.zshrc or ~/.bashrc ##############
    export PATH=/home/Aaron/mpich4/bin:$PATH
    #########################################################

The procedure to install netcdf-cxx
-------------------------------------

Install hdf5
~~~~~~~~~~~~~

.. code-block:: shell

    wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.14.tar.gz
    tar xf hdf5-1.8.14.tar.gz
    cd hdf5-1.8.14
    # ./configure --prefix=$HOME/local --enable-fortran --enable-cxx --enable-parallel --enable-unsupported (error for parallel)
    ./configure --prefix=$HOME/local --enable-fortran --enable-cxx
    make -j 32
    make check
    make install

Install netcdf-c
~~~~~~~~~~~~~~~~~

.. code-block:: shell

    wget https://github.com/Unidata/netcdf-c/archive/v4.3.3.1.tar.gz
    tar xf v4.3.3.1.tar.gz
    cd netcdf-c-4.3.3.1
    ./configure --prefix=$HOME/local --enable-netcdf-4
    make
    make check
    make install

Install netcdf-cxx
~~~~~~~~~~~~~~~~~~~~

.. code-block:: shell

    wget https://github.com/Unidata/netcdf-cxx4/archive/v4.2.1.tar.gz
    tar xf v4.2.1.tar.gz
    cd netcdf-cxx4-4.2.1
    ./configure --prefix=$HOME/local
    make
    make check
    make install

Install petsc
----------------

.. code-block:: shell

    wget ~~~~
    cd ~/petsc
    ./configure COPTFLAGS="-g -O3" --prefix=${petsc_prefix} --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpifort --download-f2cblaslapack=1 --download-triangle=1

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

