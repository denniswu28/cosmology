#!/bin/bash -l                     
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 6:00:00
#SBATCH --account=cosmology
#SBATCH --error=err/shear_calc.err 
#SBATCH --mem=100G
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

# python /hpc/group/cosmology/denniswu/cosmosis_files/shear_plotting/plot_21.py
# python /hpc/group/cosmology/denniswu/scripts/get_n_nz.py
# python /hpc/group/cosmology/denniswu/scripts/plot_covmat_default.py
python /hpc/group/cosmology/denniswu/scripts/shear_calc.py
# python /hpc/group/cosmology/denniswu/scripts/write_fits.py
