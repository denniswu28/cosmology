srun --partition cosmology --pty bash -i
source /hpc/group/cosmology/phy-lsst/cl562/lsst-pipeline/v21.0.0/loadLSST.bash;setup lsst_distrib
jupyter notebook --no-browser --port=7777 --ip=$(hostname -s)
