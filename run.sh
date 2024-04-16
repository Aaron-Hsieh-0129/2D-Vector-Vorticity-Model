#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 3-00:00
#SBATCH --mem=2000
#SBATCH -o myjob.o
#SBATCH -e myjob.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b08209006@ntu.edu.tw

rm -rf build
mkdir build
cd build/ && cmake ../ && make -j 4 && mpirun -n 1 ./vvm2d -ksp_type gmres
