#!/bin/bash
#SBATCH --account=MST113255
#SBATCH --job-name=2DVVM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=log/0514/job-%j.out
#SBATCH --error=log/0514/job-%j.err
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=zz85721@gmail.com

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
# cd build/ && cmake .. && make -j 8 && mpirun -np 1 -mca btl_base_warn_component_unused 0 -np 1 -x CUDA_VISIBLE_DEVICES=0 ./vvm2d
