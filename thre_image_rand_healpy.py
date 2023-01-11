"""
This script is to compute the 2-point count-shear correlation of the truth catalog using treecorr, with bright objects as the lenses. 
The lines reading in the catalog is adapted from Prof. Troxel's script. Ideas on generating random catalog referenced Prof. Troxel's script on 
https://github.com/matroxel/destest/blob/master/src/catalog.py#L1435

Author: Dennis Wu
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
import numpy.random as rand
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import healpy as hp
import healsparse as hsp
import skyproj
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

fr = ['Y106','J129','H158','F184']
thre_mag = 22

# Read in all objects
dc2_truth = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_tru.fits.gz')[-1].read()

# Select objects based on threshold
dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < thre_mag]

# Read in shear data
dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()

# Create new array with only shear data of bright objects
dc2_thre_shear = np.copy(dc2_truth_shear[dc2_truth_thre['ind']])

ra = dc2_truth_shear['ra']
dec = dc2_truth_shear['dec']

ra_thre = dc2_thre_shear['ra']
dec_thre = dc2_thre_shear['dec']

s1 = np.zeros(len(dc2_truth_shear))
s2 = np.zeros(len(dc2_truth_shear))

with open('s1_NG.csv') as fs1:
	reader = csv.reader(fs1, delimiter=',')
	for row in reader:
		s1 = row

with open('s2_NG.csv') as fs2:
	reader = csv.reader(fs2, delimiter=',')
	for row in reader:
		s2 = row

print("total number: ", ra.size)
print("thresholded number: ", ra_thre.size)

n_len = len(dc2_truth_thre)
del(dc2_truth)
del(dc2_thre_shear)
del(dc2_truth_thre)

# Generate a set of random ra and dec by healpy and healsprase (part of the script adpated from HealSparse tutorial)
if(rank == 0):
	nside_coverage = 32
	nside_sparse = 512
	hsp_map = hsp.HealSparseMap.make_empty(nside_coverage, nside_sparse, dtype=np.float64, sentinel=0)
	pix = hp.pixelfunc.ang2pix(nside_sparse, dc2_truth_shear['dec'], dc2_truth_shear['ra'], nest=True)
	pix_unq, counts = np.unique(pix, return_counts=True)
	hsp_map.update_values_pix(pix_unq, counts)
	fig, ax = plt.subplots()
	m = skyproj.McBrydeSkyproj(ax=ax, autorescale=False)
	_ = m.draw_hspmap(hsp_map, cmap='viridis')
	plt.savefig('hsp_map.pdf')
	plt.clf()
	hsp_map.write('hsp_map.hs', clobber=False)
else:
	hsp_map = hsp.HealSparseMap.read('hsp_map.hs')


ra_rand, dec_rand = hsp.make_uniform_randoms_fast(hsp_map, n_len)

with open('ra_rand_'+str(rank)+'_hsp.csv', 'w') as ra_rand_w:
	writer = csv.writer(ra_rand_w)
	writer.writerow(ra_rand)

with open('dec_rand_'+str(rank)+'_hsp.csv', 'w') as dec_rand_w:
	writer = csv.writer(dec_rand_w)
	writer.writerow(dec_rand)



