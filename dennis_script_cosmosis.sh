#!/bin/bash -l                     
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --account=cosmology
#SBATCH --error=err/cosmosis_0.err 
#SBATCH --mem=200G
#SBATCH -p common
#SBATCH --mail-type=all
#SBATCH --mail-user=tianrui.wu@duke.edu

cd /hpc/group/cosmology
# source loadEnv23.sh
mamba activate /hpc/group/cosmology/cosmosis_env
source cosmosis-configure

# #OpenMP settings:
export OMP_NUM_THREADS=1
# #export OMP_PLACES=threads
# #export OMP_PROC_BIND=spread
# export I_MPI_SPIN_COUNT=1
# export I_MPI_SHM_EAGER_THRESHOLD=4096
# export SLURM_CPU_BIND=none
# export SLURM_WHOLE=1

# mpirun -n 100 cosmosis --mpi /hpc/group/cosmology/denniswu/cosmosis_files/script/params_real.ini
# mpirun -n 100 cosmosis --mpi /hpc/group/cosmology/denniswu/cosmosis_files/script/params_real_2.ini
# mpirun -n 100 cosmosis --mpi /hpc/group/cosmology/denniswu/cosmosis_files/script/params_real_3.ini

# cosmosis /hpc/group/cosmology/denniswu/cosmosis_files/script/params_pseudo.ini
# cosmosis-postprocess denniswu/cosmosis_files/output/truth_truth_3x2pt/output.txt -o denniswu/cosmosis_files/output/truth_truth_3x2pt/plot -p truth_truth --derive denniswu/cosmosis_files/script/derived.py
# cosmosis-postprocess denniswu/cosmosis_files/output/roman_roman_gw/output.txt -o denniswu/cosmosis_files/output/roman_roman_gw/plot -p roman_roman --derive denniswu/cosmosis_files/script/derived.py
# cosmosis-postprocess denniswu/cosmosis_files/output/truth_truth_wl/output.txt denniswu/cosmosis_files/output/rubin_rubin_wl/output.txt denniswu/cosmosis_files/output/roman_roman_wl/output.txt -o denniswu/cosmosis_files/output/wl_3overlay/plot -p cosmic_shear --legend="truth|rubin|roman"
cosmosis-postprocess denniswu/cosmosis_files/output/truth_roman_3x2pt/output.txt -o denniswu/cosmosis_files/output/truth_roman_3x2pt/plot -p truth_roman --derive denniswu/cosmosis_files/script/derived.py
cosmosis-postprocess denniswu/cosmosis_files/output/truth_roman_wl/output.txt -o denniswu/cosmosis_files/output/truth_roman_wl/plot -p truth_roman --derive denniswu/cosmosis_files/script/derived.py
cosmosis-postprocess denniswu/cosmosis_files/output/truth_roman_gw/output.txt -o denniswu/cosmosis_files/output/truth_roman_gw/plot -p truth_roman --derive denniswu/cosmosis_files/script/derived.py
cosmosis-postprocess denniswu/cosmosis_files/output/truth_rubin_3x2pt/output.txt -o denniswu/cosmosis_files/output/truth_rubin_3x2pt/plot -p truth_rubin --derive denniswu/cosmosis_files/script/derived.py
cosmosis-postprocess denniswu/cosmosis_files/output/truth_rubin_wl/output.txt -o denniswu/cosmosis_files/output/truth_rubin_wl/plot -p truth_rubin --derive denniswu/cosmosis_files/script/derived.py
cosmosis-postprocess denniswu/cosmosis_files/output/truth_rubin_gw/output.txt -o denniswu/cosmosis_files/output/truth_rubin_gw/plot -p truth_rubin --derive denniswu/cosmosis_files/script/derived.py


# cosmosis-postprocess --help