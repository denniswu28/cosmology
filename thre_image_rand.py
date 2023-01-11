"""
This script is to compute the 2-point count-shear correlation of the truth catalog using treecorr, with bright objects as the lenses. 
The lines reading in the catalog is adapted from Prof. Troxel's script. Ideas on generating random catalog referenced Prof. Troxel's script on 
https://github.com/matroxel/destest/blob/master/src/catalog.py#L1435
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
import numpy.random as rand
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

fr = ['Y106','J129','H158','F184']

# Read in all objects
dc2_truth = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_tru.fits.gz')[-1].read(columns=['mag_J129','mag_F184','mag_H158','mag_Y106', 'ind'])

# Select objects based on threshold and thus determine the size of the random sample (10x)
dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < 23]
dc2_leftover = np.copy(dc2_truth['ind'])
del(dc2_truth)
n_rand = int(len(dc2_truth_thre)*15) # size of random sample
del(dc2_truth_thre)

# Read in shear data
dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read(columns=['ra','dec'])

# Generate a set of random ra and dec by setting subtiles
n_tile_dec = 24
n_tile_ra = 32
n_per_tile = int(n_rand/(n_tile_ra*n_tile_dec))+1
print('nr_per_tile =', n_per_tile)

# update n_rand
n_rand = n_per_tile*n_tile_ra*n_tile_dec
print('nrand =', n_rand)

# create empty ra/dec and indices array to store the final random sample
ra_rand = np.zeros(n_rand)
dec_rand = np.zeros(n_rand)
rand_ind = []
tile_cnt = np.zeros((n_tile_ra,n_tile_dec))

# set tile bounds based on min-max bounds and tile numbers
z_bound = np.linspace(-0.6156614, -0.6691307, n_tile_dec+1)
dec_bound = np.arcsin(z_bound)
ra_bound = np.linspace(0.890117, 0.977386, n_tile_ra+1)

# start uniform random sampling
while len(rand_ind) < n_rand:
	leftover_size = len(dc2_leftover)
	ind_rand = rand.randint(0,leftover_size, n_per_tile)
	for index in ind_rand:
		i = np.digitize(dc2_truth_shear[dc2_leftover[index]]['ra'], ra_bound)-1
		if(i > n_tile_ra-1 or i < 0):  # discard obj exceeding the bounds
			continue
		j = np.digitize(dc2_truth_shear[dc2_leftover[index]]['dec'], dec_bound)-1
		if(j > n_tile_dec-1 or j < 0): # discard obj exceeding the bounds
			continue
		if(tile_cnt[i, j] < n_per_tile): # take obj if correspding tile not full
			rand_ind.append(dc2_leftover[index])
			tile_cnt[i,j] += 1

# discard redundant obj
rand_ind= np.unique(rand_ind)
print("final rand catalog size:", len(rand_ind))

# using indices to take all ra/dec of selected objs
ra_rand = np.copy(dc2_truth_shear[rand_ind]['ra'])
dec_rand = np.copy(dc2_truth_shear[rand_ind]['dec'])

# output to csv file
with open('ra_rand'+'_man.csv', 'w') as ra_rand_w:
	writer = csv.writer(ra_rand_w)
	writer.writerow(ra_rand)

with open('dec_rand'+'_man.csv', 'w') as dec_rand_w:
	writer = csv.writer(dec_rand_w)
	writer.writerow(dec_rand)
