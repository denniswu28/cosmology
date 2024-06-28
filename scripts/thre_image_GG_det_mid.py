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
# dec_min = -42
# dec_max = -38
# ra_min = 51
# ra_max = 56
dec_min = -41
dec_max = -39
ra_min = 52.5
ra_max = 54.5

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# Read in and cut dectection objects
dc2_det = util.read_det_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz')
dc2_det = dc2_det[dc2_det['number']>0] # select valid
dc2_det = dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5
dc2_det = dc2_det[(dc2_det['alphawin_j2000']>ra_min)&(dc2_det['alphawin_j2000']<ra_max)&(dc2_det['deltawin_j2000']>dec_min)&(dc2_det['deltawin_j2000']<dec_max)]

# Match both catalogs
dc2_det_match, dc2_truth_match = util.get_match(dc2_det, dc2_truth)

# Extract position data
ra = dc2_det_match['alphawin_j2000']
dec = dc2_det_match['deltawin_j2000']

# Extract shear data
s1 = np.zeros(len(dc2_det_match))
s2 = np.zeros(len(dc2_det_match))

with open('data/s1_det_match_mid.csv') as fs1:
	reader = csv.reader(fs1, delimiter=',')
	for row in reader:
		s1 = row

with open('data/s2_det_match_mid.csv') as fs2:
	reader = csv.reader(fs2, delimiter=',')
	for row in reader:
		s2 = row


for i in range(20,26):
	mag_threshold = i

	# Select objects based on threshold
	dc2_det_thre = np.copy(dc2_det_match[dc2_det_match['mag_auto_'+fr[2]] < mag_threshold])
	dc2_truth_thre = np.copy(dc2_truth_match[dc2_det_match['mag_auto_'+fr[2]] < mag_threshold])

	# Extract thresholded position data
	ra_thre = dc2_det_thre['alphawin_j2000']
	dec_thre = dc2_det_thre['deltawin_j2000']

	# Extract thresholded shear data
	s1_thre = np.zeros(len(dc2_det_thre))
	s2_thre = np.zeros(len(dc2_det_thre))

	with open('data/s1_det_match_thre_mid_'+ str(mag_threshold) +'.csv') as fst1:
		reader = csv.reader(fst1, delimiter=',')
		for row in reader:
			s1_thre = row

	with open('data/s2_det_match_thre_mid_'+ str(mag_threshold) +'.csv') as fst2:
		reader = csv.reader(fst2, delimiter=',')
		for row in reader:
			s2_thre = row


	# Generate the treecorr catalogs
	k = 10
	cat = treecorr.Catalog(ra = ra, dec = dec, g1 = s1, g2 = s2, ra_units='degrees', dec_units='degrees', npatch = k)
	cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, g1 = s1_thre, g2 = s2_thre, ra_units='degrees', dec_units='degrees', patch_centers=cat.patch_centers)

	# Construct the count-shear corr object
	gg = treecorr.GGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='marked_bootstrap')
	gg_thre = treecorr.GGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='marked_bootstrap')
	gg_cross = treecorr.GGCorrelation(min_sep=1, max_sep=100, nbins=20, sep_units='arcmin', var_method='marked_bootstrap')

	# Compute the shear-shear corr and compensate with rand background
	gg.process(cat)
	gg_thre.process(cat_thre)
	gg_cross.process(cat_thre, cat)
	# gg.write("data/gg_corr_det_all.fits")
	# gg_thre.write("data/gg_corr_det_"+ str(mag_threshold) +".fits")
	# gg_cross.write("data/gg_cross_det_"+ str(mag_threshold) +".fits")
	gg.write("data/gg_corr_det_mid_all.fits")
	gg_thre.write("data/gg_corr_det_mid_"+ str(mag_threshold) +".fits")
	gg_cross.write("data/gg_cross_det_mid_"+ str(mag_threshold) +".fits")




