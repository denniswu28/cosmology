"""
This script is to compute the 2-point point-shear correlation of the threshholded truth catalog using treecorr. The lines reading in the catalog is adapted from Prof. Troxel's script.

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
# from mpi4py import MPI

# self.comm = MPI.COMM_WORLD
# self.rank = self.comm.Get_rank()
# self.size = self.comm.Get_size()

fr = ['Y106','J129','H158','F184']
dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# Read in shear data
dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()
# dc2_truth_shear = dc2_truth_shear[np.logical_and(np.logical_and(dc2_truth_shear['ra'] < np.deg2rad(ra_max), dc2_truth_shear['ra'] > np.deg2rad(ra_min)),np.logical_and(dc2_truth_shear['dec'] < np.deg2rad(dec_max), dc2_truth_shear['dec'] > np.deg2rad(dec_min)))]
dc2_truth_shear_limited = np.copy(dc2_truth_shear[dc2_truth['ind']])

# Extract position data
ra = dc2_truth_shear_limited['ra']
dec = dc2_truth_shear_limited['dec']

# Extract shear data
s1 = np.zeros(len(dc2_truth_shear_limited))
s2 = np.zeros(len(dc2_truth_shear_limited))

del(dc2_truth_shear_limited)

with open('data/s1.csv') as fs1:
	reader = csv.reader(fs1, delimiter=',')
	for row in reader:
		s1 = row

with open('data/s2.csv') as fs2:
	reader = csv.reader(fs2, delimiter=',')
	for row in reader:
		s2 = row


for i in range(20,26):
	mag_threshold = i

	# Select objects based on threshold
	dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < mag_threshold]

	# Extract position and shear data
	ra_thre = dc2_truth_thre['ra']
	dec_thre = dc2_truth_thre['dec']

	del(dc2_truth_thre)

	# Read the random catalog to arrays
	with open('data/ra_rand_gal_'+ str(mag_threshold) +'.csv', newline='') as rafile:
		ra_rand = np.float_(list(csv.reader(rafile, delimiter=',')))
	ra_rand = ra_rand[0]

	with open('data/dec_rand_gal_'+ str(mag_threshold) +'.csv', newline='') as decfile:
		dec_rand = np.float_(list(csv.reader(decfile, delimiter=',')))
	dec_rand = dec_rand[0]

	# Generate the treecorr catalogs
	k = 10
	cat = treecorr.Catalog(ra = ra, dec = dec, g1 = s1, g2 = s2, ra_units='radians', dec_units='radians', npatch = k)
	cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, ra_units='degrees', dec_units='degrees', patch_centers=cat.patch_centers)
	cat_rand = treecorr.Catalog(ra = ra_rand, dec = dec_rand, ra_units='degrees', dec_units='degrees', patch_centers=cat.patch_centers)

	# Construct the count-shear corr object
	ng = treecorr.NGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
	ng_rand = treecorr.NGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')

	# Compute the count-shear corr and compensate with rand background
	ng.process(cat_thre, cat)
	ng_rand.process(cat_rand, cat)
	ng.calculateXi(rg = ng_rand)
	ng.write("data/ng_corr_"+ str(mag_threshold) +".fits")
	ng_rand.write("data/ng_rand_"+ str(mag_threshold) +".fits")

	# # plot the xi functions
	# r = np.exp(ng.meanlogr)
	# xi = ng.xi
	# sig = np.sqrt(ng.varxi)

	# plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='blue', label = r'$\xi_(\theta)$')

	# plt.xscale('log')
	# plt.yscale('log')
	# #plt.xlim([0.01,250])

	# plt.xlabel(r'$\theta$ (arcmin)')
	# plt.ylabel(r'$\xi$')
	# plt.legend()
	# plt.title("2pt correlation of galaxy shear in CosmoDC2")

	# # Save both functions in a pdf file
	# plt.savefig('pdf/corr_func_threshold_NG.pdf')
	# plt.clf()























