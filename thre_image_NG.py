"""
This script is to compute the 2-point shear-shear correlation of the threshholded truth catalog using treecorr. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import treecorr


fr = ['Y106','J129','H158','F184']

# Read in all objects
dc2_truth = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_tru.fits.gz')[-1].read()

# Select objects based on threshold
dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < 23]

# Read in shear data
dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()

# Create new array with only data of bright objects
dc2_thre_shear = np.copy(dc2_truth_shear[dc2_truth_thre['ind']])

# Extract position and shear data
ra = dc2_truth_shear['ra']
dec = dc2_truth_shear['dec']

ra_thre = dc2_thre_shear['ra']
dec_thre = dc2_thre_shear['dec']

s1 = np.zeros(len(dc2_truth_shear))
s2 = np.zeros(len(dc2_truth_shear))

g1 = dc2_truth_shear['g1']
g2 = dc2_truth_shear['g2']


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

del(dc2_truth)
del(dc2_thre_shear)
del(dc2_truth_thre)

# Read the random catalog to arrays
ra_rand = np.array([])
dec_rand = np.array([])

with open('ra_rand_man.csv', newline='') as rafile:
    ra_rand = list(csv.reader(rafile, delimiter=','))

with open('dec_rand_man.csv', newline='') as decfile:
    dec_rand = list(csv.reader(decfile, delimiter=','))

# Generate the treecorr catalogs
cat = treecorr.Catalog(ra = ra, dec = dec, g1 = s1, g2 = s2, ra_units='radians', dec_units='radians', npatch = 20)
cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, ra_units='radians', dec_units='radians', npatch = 20)
cat_rand = treecorr.Catalog(ra = ra_rand, dec = dec_rand, ra_units='radians', dec_units='radians', npatch = 20)

# Construct the count-shear corr object
ng = treecorr.NGCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin', var_method='jackknife')
ng_rand = treecorr.NGCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin', var_method='jackknife')

# Compute the count-shear corr and compensate with rand background
ng.process(cat_thre, cat)
ng_rand.process(cat_rand, cat)
ng.calculateXi(rg = ng_rand)
ng.write("ng_corr_reduced.fits")

# plot the xi functions
r = np.exp(ng.meanlogr)
xi = ng.xi
sig = np.sqrt(ng.varxi)

plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='blue', label = r'$\xi_(\theta)$')

plt.xscale('log')
plt.yscale('log')
#plt.xlim([0.01,250])

plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\xi$')
plt.legend()
plt.title("2-point count-shear correlation")

# Save both functions in a pdf file
plt.savefig('corr_func_threshold.pdf')
plt.clf()


# plot the xi functions
r = np.exp(ng_rand.meanlogr)
xi = ng_rand.xi
sig = np.sqrt(ng_rand.varxi)

plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='blue', label = r'$\xi_(\theta)$')

#plt.xscale('log')
#plt.yscale('log')
#plt.xlim([0.01,250])

plt.xlabel(r'$\theta$ (arcmin)')
plt.ylabel(r'$\xi$')
plt.legend()
plt.title("2-point count-shear correlation")

# Save both functions in a pdf file
plt.savefig('corr_func_threshold_rand.pdf')
plt.clf()




















