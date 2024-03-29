#!/bin/bash -l                     
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -t 20:00:00
#SBATCH --account=cosmology
#SBATCH --error=err/cosmosis.err 
#SBATCH --mem=100G
#SBATCH -p cosmology
#SBATCH --mail-type=all
#SBATCH --mail-user=tianrui.wu@duke.edu

cd /hpc/group/cosmology
source /hpc/group/cosmology/loadEnv23.sh
# source /hpc/group/cosmology/loadEnv.sh
source /hpc/group/cosmology/cosmosis-configure

#OpenMP settings:
export OMP_NUM_THREADS=1
#export OMP_PLACES=threads
#export OMP_PROC_BIND=spread
export I_MPI_SPIN_COUNT=1
export I_MPI_SHM_EAGER_THRESHOLD=4096
export SLURM_CPU_BIND=none
export SLURM_WHOLE=1

# cosmosis /hpc/group/cosmology/denniswu/cosmosis_files/script/params.ini
cosmosis-postprocess denniswu/cosmosis_files/output/CosmoDC2_shear_01_31/output.txt denniswu/cosmosis_files/output/CosmoDC2/output.txt -o denniswu/cosmosis_files/output/CosmoDC2_shear_02_07/plot -p CosmoDC2_shear_02_07 --legend="det|truth"
# cosmosis-postprocess --help