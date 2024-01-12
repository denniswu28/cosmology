"""
This script is to generate a randomized, uniform catalog from the truth catalog that covers a rectangular survey region.  
The lines reading in the catalog is adapted from Prof. Troxel's script. Ideas on generating random catalog referenced Prof. Troxel's script on 
https://github.com/matroxel/destest/blob/master/src/catalog.py#L1435
"""

# All the libs needed
import csv
import numpy as np
from numpy.random import default_rng
import utilities as util
rng = default_rng()

fr = ['Y106','J129','H158','F184']
mag_threshold = 25
dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# Select objects based on threshold and thus determine the size of the random sample (15x)
dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < mag_threshold]
n_rand = int(len(dc2_truth_thre)*15) # size of random sample
del(dc2_truth_thre)

# Generate a set of random ra and dec by setting subtiles
n_tile_dec = 9
n_tile_ra = 12
n_per_tile = int(n_rand/(n_tile_ra*n_tile_dec))+1
print('n_per_tile =', n_per_tile)

# update n_rand
n_rand = n_per_tile*n_tile_ra*n_tile_dec
print('n_rand =', n_rand)

# create empty ra/dec and indices array to store the final random sample
ra_rand = []
dec_rand = []
rand_ind = []
tile_cnt = np.zeros((n_tile_ra,n_tile_dec))

# set tile bounds based on min-max bounds and tile numbers
z_min = np.sin(np.deg2rad(dec_min))
z_max = np.sin(np.deg2rad(dec_max))
z_bound = np.linspace(z_min, z_max, n_tile_dec+1)
dec_bound = np.rad2deg(np.arcsin(z_bound))
ra_bound = np.linspace(ra_min, ra_max, n_tile_ra+1)
ind_size = len(dc2_truth)

# start uniform random sampling
while len(rand_ind) < n_rand:
	ind_rand = rng.integers(0,ind_size, n_per_tile)
	for index in ind_rand:
		i = np.digitize(dc2_truth[index]['ra'], ra_bound)-1
		print("i =", i, "; ra = ", dc2_truth[index]['ra'])
		if(i > n_tile_ra-1 or i < 0):  # discard obj exceeding the bounds
			continue
		j = np.digitize(dc2_truth[index]['dec'], dec_bound)-1
		print("j =", j, "; dec = ", dc2_truth[index]['dec'])
		if(j > n_tile_dec-1 or j < 0): # discard obj exceeding the bounds
			continue
		if(tile_cnt[i, j] < n_per_tile): # take obj if correspding tile not full
			rand_ind.append(index)
			tile_cnt[i,j] += 1

# discard redundant obj
rand_ind = np.unique(rand_ind)
print("final rand catalog size:", len(rand_ind))

# using indices to take all ra/dec of selected objs
ra_rand = np.copy(dc2_truth[rand_ind]['ra'])
dec_rand = np.copy(dc2_truth[rand_ind]['dec'])

# output to csv file
with open('data/ra_rand_gal_'+ str(mag_threshold) + '.csv', 'w') as ra_rand_w:
	writer = csv.writer(ra_rand_w)
	writer.writerow(ra_rand)

with open('data/dec_rand_gal_'+ str(mag_threshold) + '.csv', 'w') as dec_rand_w:
	writer = csv.writer(dec_rand_w)
	writer.writerow(dec_rand)