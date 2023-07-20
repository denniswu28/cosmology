#!/bin/bash -l                     
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -t 1:30:00
#SBATCH --account=cosmology
#SBATCH --error=err/write_fits.err 
#SBATCH --mem=200G
#SBATCH -p cosmology
#SBATCH --mail-type=all
#SBATCH --mail-user=tianrui.wu@duke.edu

source /hpc/group/cosmology/loadEnv.sh

#OpenMP settings:
export OMP_NUM_THREADS=1
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread
export I_MPI_SPIN_COUNT=1
export I_MPI_SHM_EAGER_THRESHOLD=4096
export SLURM_CPU_BIND=none
export SLURM_WHOLE=1

python /hpc/group/cosmology/denniswu/write_fits.py

