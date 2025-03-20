#!/bin/bash
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --nodelist=mogamd
#SBATCH -o log/ShearNoCouple.o
#SBATCH -e log/ShearNoCouple.e
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b08209006@ntu.edu.tw

rm -rf build
mkdir build

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads
echo $OMP_NUM_THREADS

# cd build/ && cmake ../ && make -j 4 && mpirun -n 1 ./vvm2d -ksp_type gmres
cd build/ && cmake ../ && make -j 8 && ./vvm2d
