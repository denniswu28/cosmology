"""
This script is to compute the 2-point count-count correlation of the threshholded truth catalog using treecorr with Landy-Szalay estimator. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import csv
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import treecorr
import utilities as util

fr = ['Y106','J129','H158','F184']
mag_threshold = 23
dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]


# Select objects based on threshold
dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < mag_threshold]

# Extract position and shear data
ra_thre = dc2_truth_thre['ra']
dec_thre = dc2_truth_thre['dec']

# Release RAM
del(dc2_truth)
del(dc2_truth_thre)

# Read the random catalog to arrays
with open('ra_rand_gal.csv', newline='') as rafile:
    ra_rand = np.float_(list(csv.reader(rafile, delimiter=',')))
ra_rand = ra_rand[0]

with open('dec_rand_gal.csv', newline='') as decfile:
    dec_rand = np.float_(list(csv.reader(decfile, delimiter=',')))
dec_rand = dec_rand[0]

# Generate the treecorr catalogs and patches with k-means
k = 10
cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, ra_units='degrees', dec_units='degrees', npatch = k)
cat_rand = treecorr.Catalog(ra = ra_rand, dec = dec_rand, ra_units='degrees', dec_units='degrees', patch_centers=cat_thre.patch_centers)

# Construct the count-shear corr object
dd = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
dr = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin')
rr = treecorr.NNCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin')

# Compute the count-count corr
dd.process(cat_thre)
dr.process(cat_thre, cat_rand)
rr.process(cat_rand)
dd.calculateXi(rr=rr, dr = dr)
dd.write("nn_corr_match_gal.fits")

# plot the corr functions for NN
r = np.exp(dd.meanlogr)
xi = dd.xi
sig = np.sqrt(dd.varxi)

plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='blue', label = r'$\xi_(\theta)$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\xi$')
plt.legend()
plt.title("2-point NN correlation")

plt.savefig('corr_func_threshold_NN.pdf')
plt.clf()





