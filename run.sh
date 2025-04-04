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

# CPU only
cd build/ && cmake ../ && make -j 8 && mpirun -np 1 ./vvm2d

# wiht GPU
# cd build/ && cmake ../ && make -j 8 && mpirun -np 1 -mca btl_base_warn_component_unused 0 -np 1 -x CUDA_VISIBLE_DEVICES=0 ./vvm2d
