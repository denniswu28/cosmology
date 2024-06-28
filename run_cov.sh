#!/bin/bash -l                     
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 6:00:00
#SBATCH --account=cosmology
#SBATCH --error=err/shear_calc.err 
#SBATCH --mem=100G
#SBATCH -p common
#SBATCH --mail-type=all
#SBATCH --mail-user=tianrui.wu@duke.edu

source /hpc/group/cosmology/loadEnv23.sh
cd /hpc/group/cosmology/repos/CosmoCov/covs
# make covs
# ./cov 1 /hpc/group/cosmology/denniswu/covmat/cov_pseudo_roman.ini
for i in {1..1770}; do ./cov $i /hpc/group/cosmology/denniswu/covmat/cov_pseudo_roman.ini; done


