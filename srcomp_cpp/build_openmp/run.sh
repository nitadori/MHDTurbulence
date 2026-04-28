#!/bin/bash
#SBATCH -J kokkos_omp
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -c 8
#SBATCH -t 01:00:00
#SBATCH -p genoa

module load mpi/openmpi-x86_64
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# MPI rank ごとに 8 core を排他的に割り当て
#srun --cpu-bind=cores ./Simulation.x # 動かない
mpirun -np ${SLURM_NTASKS} ./Simulation.x
