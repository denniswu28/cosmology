#!/bin/bash -l                     
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --account=cosmology
#SBATCH --error=err/shear_calc.err 
#SBATCH --mem=100G
#SBATCH -p cosmology
#SBATCH --mail-type=all
#SBATCH --mail-user=tianrui.wu@duke.edu

# source /hpc/group/cosmology/loadEnv.sh
# source /hpc/group/cosmology/loadEnv23.sh
# source cosmosis-configure
cd /hpc/group/cosmology
mamba activate /hpc/group/cosmology/mambaforge/envs/simdevel
export PYTHONPATH=${PYTHONPATH}:/hpc/group/cosmology/repos/2point/
# source cosmosis-configure



# #OpenMP settings:
# export OMP_NUM_THREADS=1
# #export OMP_PLACES=threads
# #export OMP_PROC_BIND=spread
# export I_MPI_SPIN_COUNT=1
# export I_MPI_SHM_EAGER_THRESHOLD=4096
# export SLURM_CPU_BIND=none
# export SLURM_WHOLE=1

# python /hpc/group/cosmology/denniswu/cosmosis_files/shear_plotting/plot_21.py
# python /hpc/group/cosmology/denniswu/scripts/get_n_nz.py
# python /hpc/group/cosmology/denniswu/scripts/plot_covmat_default.py
# python /hpc/group/cosmology/denniswu/scripts/thre_image_GG_read.py
# python /hpc/group/cosmology/denniswu/scripts/thre_image_GG_det.py
# python /hpc/group/cosmology/denniswu/scripts/thre_image_rand_match.py
# python /hpc/group/cosmology/denniswu/scripts/shear_calc.py
# python /hpc/group/cosmology/denniswu/scripts/write_fits.py
# python /hpc/group/cosmology/denniswu/scripts/test.py
# python /hpc/group/cosmology/denniswu/scolnic/quasar.py
# python /hpc/group/cosmology/denniswu/scripts/shear_calc_match.py
# python /hpc/group/cosmology/denniswu/scripts/match_plot.py
# python /hpc/group/cosmology/denniswu/scripts/write_fit_tomo_pseudo.py
python /hpc/group/cosmology/denniswu/scripts/write_fit_tomo.py

# cd /hpc/group/cosmology/repos/CosmoCov/covs
# make covs

# for i in {1..1485}; do ./cov $i /hpc/group/cosmology/denniswu/covmat/cov_pseudo.ini; done