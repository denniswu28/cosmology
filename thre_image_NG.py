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
import galsim
import os
import treecorr
import healpy as hp

fr = ['Y106','J129','H158','F184']


# Read in all objects
dc2_truth = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_tru.fits.gz')[-1].read()

# Select objects based on threshold
dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < 22]

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

# Generate a set of random ra and dec by setting subtiles

n_rand = len(dc2_truth_thre)*10
n_tile_dec = 3
n_tile_ra = 4
n_per_tile = int(n_rand/(n_tile_ra*n_tile_dec))+1
print('nr_per_tile =', n_per_tile)
n_rand = n_per_tile*n_tile_ra*n_tile_dec
print('nrand =', n_rand)

ra_rand = np.zeros(n_rand)
dec_rand = np.zeros(n_rand)

z_bound = np.linspace(-0.615661475, -0.669130606, n_tile_dec+1)
dec_bound = np.arcsin(z_bound)
ra_bound = np.linspace(0.890117, 0.977385, n_tile_ra+1)

rand_ind = []
num = 0
for i in range(n_tile_ra):
	for j in range(n_tile_dec):
		num = num + n_per_tile
		while len(rand_ind) < num:
			ind_rand = rand.choice(dc2_truth['ind'], size = n_per_tile, replace=False)
			for index in ind_rand:
				if(dc2_truth_shear[index]['ra'] <= ra_bound[i+1] and dc2_truth_shear[index]['ra'] >= ra_bound[i]):
					if(dc2_truth_shear[index]['dec'] <= dec_bound[j+1] and dc2_truth_shear[index]['dec'] >= dec_bound[j]):
						rand_ind.append(index)
		rand_ind = rand_ind[:num]
		
dc2_rand_shear = np.copy(dc2_truth_shear[rand_ind])
ra_rand = dc2_rand_shear['ra']
dec_rand = dc2_rand_shear['dec']

# Generate the treecorr catalog
cat = treecorr.Catalog(ra = ra, dec = dec, g1 = s1, g2 = s2, ra_units='radians', dec_units='radians')
cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, ra_units='radians', dec_units='radians')
cat_rand = treecorr.Catalog(ra = ra_rand, dec = dec_rand, ra_units='radians', dec_units='radians')
cat.write('cat.dat')
cat_thre.write('cat_thre.dat')
cat_rand.write('cat_rand.dat')

# Construct the count-shear corr object
ng = treecorr.NGCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin')

# Compute the count-shear corr
ng.process(cat_thre, cat)
ng.write("ng_corr.fits")
ng.calculateXi(rg = cat_rand)
ng.write("ng_corr_reduced.fits")

# Plot the xi functions
r = np.exp(ng.meanlogr)
xi = ng.xi
sig = np.sqrt(ng.varxi)

#plt.plot(r, xi, color='blue', linestyle="",marker=".", label = r'$\xi_(\theta)$')
plt.errorbar(r, xi, yerr=sig, color='blue', linestyle="",marker=".", label = r'$\xi_(\theta)$')


plt.xscale('log')
plt.yscale('log')
plt.xlim([0,250])
plt.xlabel(r'$\theta$ (arcmin)')

plt.legend()
plt.ylabel(r'$\xi$')

# Save both functions in a pdf file
plt.savefig('corr_func_threshold.pdf')
plt.clf()





#dc2_truth_ra_bound = dc2_truth_shear[dc2_truth_shear['ra'] >= ra_bound[i]] 
#print("max ra:", np.max(dc2_truth_shear['ra']))
#print("max dec:", np.max(dc2_truth_shear['dec']))
#dc2_truth_ra_bound = dc2_truth_ra_bound[dc2_truth_ra_bound['ra'] <= ra_bound[i+1]]
#dc2_truth_bound = dc2_truth_ra_bound[dc2_truth_ra_bound['dec'] >= dec_bound[j]]
#dc2_truth_bound = dc2_truth_bound[dc2_truth_bound['dec'] <= dec_bound[j+1]]
#print(dc2_truth_bound.shape)
#print("number: ", len(dc2_truth_bound))
#ind_rand = np.random.choice(dc2_truth_bound['ind'], n_per_tile, replace=False)
#for k in range(n_per_tile):
#dc2_truth_rand = np.copy(dc2_truth_shear[ind_rand])
#ra_rand[n_tile_dec*n_per_tile*i+n_per_tile*j] = dc2_truth_rand[k, 'ra']
#dec_rand[n_tile_dec*n_per_tile*i+n_per_tile*j] = dc2_truth_rand[k, 'dec']














