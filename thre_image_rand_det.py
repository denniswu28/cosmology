"""
This script is to compute the 2-point count-shear correlation of the detection catalog using treecorr, with bright objects as the lenses. 
The lines reading in the catalog is adapted from Prof. Troxel's script. Ideas on generating random catalog referenced Prof. Troxel's script on 
https://github.com/matroxel/destest/blob/master/src/catalog.py#L1435
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
from numpy.random import default_rng
import glob
rng = default_rng()

fr = ['Y106','J129','H158','F184']
mag_threshold = 22

# Read in all dectection objects
# dc2_det = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_det.fits.gz')[-1].read(columns=['mag_auto_J129','mag_auto_F184','mag_auto_H158','mag_auto_Y106', 'number', 'alphawin_j2000', 'deltawin_j2000','flux_auto','fluxerr_auto','flags'])

start  = 0
dc2_det = None
for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz'))):
	try:
		tmp = fio.FITS(f)[-1].read(columns=['mag_auto_J129','mag_auto_F184','mag_auto_H158','mag_auto_Y106', 'number', 'alphawin_j2000', 'deltawin_j2000','flux_auto','fluxerr_auto','flags'])
		if dc2_det is None:
			dc2_det = np.zeros(100000000,dtype=tmp.dtype)
		for col in dc2_det.dtype.names:
			dc2_det[col][start:start+len(tmp)] = tmp[col]
		start+=len(tmp)
	except:
		print('-----fail'+f)
		pass



dc2_det=dc2_det[dc2_det['number']>0] # select positive
dc2_det=dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5
dc2_det = dc2_det[ (dc2_det['alphawin_j2000']>51)&(dc2_det['alphawin_j2000']<56)&(dc2_det['deltawin_j2000']>-42)&(dc2_det['deltawin_j2000']<-38)]

# Select objects based on threshold and thus determine the size of the random sample (15x)
dc2_det_thre = dc2_det[dc2_det['mag_auto_'+fr[2]] < mag_threshold]
n_rand = int(len(dc2_det_thre)*15) # size of random sample
del(dc2_det_thre)

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
dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56
z_min = np.sin(np.deg2rad(dec_min))
z_max = np.sin(np.deg2rad(dec_max))

z_bound = np.linspace(z_min, z_max, n_tile_dec+1)
dec_bound = np.rad2deg(np.arcsin(z_bound))
ra_bound = np.linspace(ra_min, ra_max, n_tile_ra+1)

ind_size = len(dc2_det)

# start uniform random sampling
# effcnt = 0
while len(rand_ind) < n_rand:
	ind_rand = rng.integers(0,ind_size, n_per_tile)
	for index in ind_rand:
		i = np.digitize(dc2_det[index]['alphawin_j2000'], ra_bound)-1
		#print("i =", i, "; ra = ", dc2_det[index]['alphawin_j2000'])
		if(i > n_tile_ra-1 or i < 0):  # discard obj exceeding the bounds
			continue
		j = np.digitize(dc2_det[index]['deltawin_j2000'], dec_bound)-1
		#print("j =", j, "; dec = ", dc2_det[index]['deltawin_j2000'])
		if(j > n_tile_dec-1 or j < 0): # discard obj exceeding the bounds
			continue
		if(tile_cnt[i, j] < n_per_tile): # take obj if correspding tile not full
			rand_ind.append(index)
			tile_cnt[i,j] += 1


# discard redundant obj
rand_ind = np.unique(rand_ind)
print("final rand catalog size:", len(rand_ind))

# using indices to take all ra/dec of selected objs
ra_rand = np.copy(dc2_det[rand_ind]['alphawin_j2000'])
dec_rand = np.copy(dc2_det[rand_ind]['deltawin_j2000'])

# output to csv file
with open('ra_rand_det_22.csv', 'w') as ra_rand_w:
	writer = csv.writer(ra_rand_w)
	writer.writerow(ra_rand)

with open('dec_rand_det_22.csv', 'w') as dec_rand_w:
	writer = csv.writer(dec_rand_w)
	writer.writerow(dec_rand)
