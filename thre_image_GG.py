"""
This script is to compute the 2-point shear-shear correlation of the threshholded truth catalog using treecorr. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
import treecorr
import utilities as util

fr = ['Y106','J129','H158','F184']

fr_num = 2
dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()
dc2_truth_shear = dc2_truth_shear[np.logical_and(np.logical_and(dc2_truth_shear['ra'] < np.deg2rad(ra_max), dc2_truth_shear['ra'] > np.deg2rad(ra_min)),np.logical_and(dc2_truth_shear['dec'] < np.deg2rad(dec_max), dc2_truth_shear['dec'] > np.deg2rad(dec_min)))]


# Extract position data
ra = dc2_truth_shear['ra']
dec = dc2_truth_shear['dec']

# Extract shear data
s1 = np.zeros(len(dc2_truth_shear))
s2 = np.zeros(len(dc2_truth_shear))

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
	dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[fr_num]] < mag_threshold]

	# Read in thresholded shear data
	dc2_thre_shear = np.copy(dc2_truth_shear[dc2_truth_thre['ind']])
	dc2_thre_shear = dc2_thre_shear[np.logical_and(np.logical_and(dc2_thre_shear['ra'] < np.deg2rad(ra_max), dc2_thre_shear['ra'] > np.deg2rad(ra_min)),np.logical_and(dc2_thre_shear['dec'] < np.deg2rad(dec_max), dc2_thre_shear['dec'] > np.deg2rad(dec_min)))]

	# Extract thresholded position data
	ra_thre = dc2_thre_shear['ra']
	dec_thre = dc2_thre_shear['dec']

	# Extract thresholded shear data
	s1_thre = np.zeros(len(dc2_thre_shear))
	s2_thre = np.zeros(len(dc2_thre_shear))

	with open('data/s1_thre_'+ str(mag_threshold) +'.csv') as fst1:
		reader = csv.reader(fst1, delimiter=',')
		for row in reader:
			s1_thre = row

	with open('data/s2_thre_'+ str(mag_threshold) +'.csv') as fst2:
		reader = csv.reader(fst2, delimiter=',')
		for row in reader:
			s2_thre = row


	# Generate the treecorr catalogs
	k = 10
	cat = treecorr.Catalog(ra = ra, dec = dec, g1 = s1, g2 = s2, ra_units='radians', dec_units='radians', npatch = k)
	cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, g1 = s1_thre, g2 = s2_thre, ra_units='radians', dec_units='radians', patch_centers=cat.patch_centers)

	# Construct the count-shear corr object
	gg = treecorr.GGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
	gg_thre = treecorr.GGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')
	gg_cross = treecorr.GGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='jackknife')

	# Compute the count-shear corr and compensate with rand background
	gg.process(cat)
	gg_thre.process(cat_thre)
	gg_cross.process(cat_thre, cat)
	gg.write("data/gg_corr_all.fits")
	gg_thre.write("data/gg_corr_"+ str(mag_threshold) +".fits")
	gg_cross.write("data/gg_cross_"+ str(mag_threshold) +".fits")













