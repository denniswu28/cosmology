#!/bin/bash -l                     
#SBATCH -N 2
#SBATCH -n 100
#SBATCH -t 12:00:00
#SBATCH --account=cosmology
#SBATCH --error=err/cosmosis.err 
#SBATCH --mem=100G
#SBATCH -p cosmology
#SBATCH --mail-type=all
#SBATCH --mail-user=tianrui.wu@duke.edu

cd /hpc/group/cosmology
# source loadEnv23.sh
mamba activate /hpc/group/cosmology/cosmosis_env
source cosmosis-configure

# #OpenMP settings:
# export OMP_NUM_THREADS=1
# #export OMP_PLACES=threads
# #export OMP_PROC_BIND=spread
# export I_MPI_SPIN_COUNT=1
# export I_MPI_SHM_EAGER_THRESHOLD=4096
# export SLURM_CPU_BIND=none
# export SLURM_WHOLE=1

mpirun -n 100 cosmosis --mpi /hpc/group/cosmology/denniswu/cosmosis_files/script/params_real.ini

# cosmosis /hpc/group/cosmology/denniswu/cosmosis_files/script/params_pseudo.ini
# cosmosis-postprocess denniswu/cosmosis_files/output/truth_truth_wl/output.txt -o denniswu/cosmosis_files/output/truth_truth_real/plot -p truth_truth --derive denniswu/cosmosis_files/script/derived.py
# cosmosis-postprocess denniswu/cosmosis_files/output/rubin_rubin_wl/output.txt -o denniswu/cosmosis_files/output/rubin_rubin_real/plot -p rubin_rubin --derive denniswu/cosmosis_files/script/derived.py
# cosmosis-postprocess denniswu/cosmosis_files/output/truth_truth_wl/output.txt denniswu/cosmosis_files/output/rubin_rubin_wl/output.txt -o denniswu/cosmosis_files/output/wl/plot -p cosmic_shear --legend="truth|rubin"
# cosmosis-postprocess --help