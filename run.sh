#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 3-00:00
#SBATCH --mem=2000
#SBATCH -o myjob.o
#SBATCH -e myjob.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b08209006@ntu.edu.tw

rm -rf build outputs
mkdir build
cd build/ && cmake ../ && make && ./vvm2d
