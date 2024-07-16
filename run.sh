#!/bin/bash
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -o myjob.o
#SBATCH -e myjob.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b08209006@ntu.edu.tw

rm -rf build
mkdir build
# cd build/ && cmake ../ && make -j 4 && mpirun -n 1 ./vvm2d -ksp_type gmres
cd build/ && cmake ../ && make -j 4 && ./vvm2d
