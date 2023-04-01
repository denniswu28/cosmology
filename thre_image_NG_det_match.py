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
import utilities as util

fr = ['Y106','J129','H158','F184']
mag_threshold = 20
dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56

# Read in and cut dectection objects
dc2_det = util.read_det_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz')
dc2_det = dc2_det[dc2_det['number']>0] # select valid
dc2_det = dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5
dc2_det = dc2_det[(dc2_det['alphawin_j2000']>ra_min)&(dc2_det['alphawin_j2000']<ra_max)&(dc2_det['deltawin_j2000']>dec_min)&(dc2_det['deltawin_j2000']<dec_max)]

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# Match both catalogs
dc2_det_match, dc2_truth_match = util.get_match(dc2_det, dc2_truth)

# Select objects based on threshold
dc2_det_thre = dc2_det_match[dc2_det_match['mag_auto_'+fr[2]] < mag_threshold]

# Read in shear data
dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()
dc2_truth_shear = dc2_truth_shear[np.logical_and(np.logical_and(dc2_truth_shear['ra'] < np.deg2rad(ra_max), dc2_truth_shear['ra'] > np.deg2rad(ra_min)),np.logical_and(dc2_truth_shear['dec'] < np.deg2rad(dec_max), dc2_truth_shear['dec'] > np.deg2rad(dec_min)))]

# Extract position and shear data
ra = dc2_truth_shear['ra']
dec = dc2_truth_shear['dec']

ra_thre = dc2_det_thre['alphawin_j2000']
dec_thre = dc2_det_thre['deltawin_j2000']

s1 = np.zeros(len(dc2_truth_shear))
s2 = np.zeros(len(dc2_truth_shear))

with open('data/s1_NG.csv') as fs1:
	reader = csv.reader(fs1, delimiter=',')
	for row in reader:
		s1 = row

with open('data/s2_NG.csv') as fs2:
	reader = csv.reader(fs2, delimiter=',')
	for row in reader:
		s2 = row

print("total number: ", ra.size)
print("thresholded number: ", ra_thre.size)

del(dc2_truth)
del(dc2_det)
del(dc2_det_thre)

# Read the random catalog to arrays
with open('data/ra_rand_det_match_'+ str(mag_threshold) +'.csv', newline='') as rafile:
    ra_rand = np.float_(list(csv.reader(rafile, delimiter=',')))
ra_rand = ra_rand[0]

with open('data/dec_rand_det_match_'+ str(mag_threshold) +'.csv', newline='') as decfile:
    dec_rand = np.float_(list(csv.reader(decfile, delimiter=',')))
dec_rand = dec_rand[0]

# Generate the treecorr catalogs
k = 10
cat = treecorr.Catalog(ra = ra, dec = dec, g1 = s1, g2 = s2, ra_units='radians', dec_units='radians', npatch = k)
cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, ra_units='degrees', dec_units='degrees', patch_centers=cat.patch_centers)
cat_rand = treecorr.Catalog(ra = ra_rand, dec = dec_rand, ra_units='degrees', dec_units='degrees', patch_centers=cat.patch_centers)

# Construct the count-shear corr object
ng = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')
ng_rand = treecorr.NGCorrelation(min_sep=0.6, max_sep=60, nbins=20, sep_units='arcmin', var_method='jackknife')

# Compute the count-shear corr and compensate with rand background
ng.process(cat_thre, cat)
ng_rand.process(cat_rand, cat)
ng.calculateXi(rg = ng_rand)
ng.write("data/ng_corr_matched_det_"+ str(mag_threshold) +".fits")

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
plt.title("2pt correlation of galaxy shear in CosmoDC2 det")

# Save both functions in a pdf file
plt.savefig('pdf/corr_func_threshold_NG_det.pdf')
plt.clf()





















