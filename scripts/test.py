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
import galsim
import pylab
import healpy as hp
import os
import treecorr
import glob
import utilities as util


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

# Read in and cut dectection objects
dc2_det = util.read_det_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz')
dc2_det = dc2_det[dc2_det['number']>0] # select valid
dc2_det = dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5
dc2_det = dc2_det[(dc2_det['alphawin_j2000']>ra_min)&(dc2_det['alphawin_j2000']<ra_max)&(dc2_det['deltawin_j2000']>dec_min)&(dc2_det['deltawin_j2000']<dec_max)]

# Match both catalogs
dc2_det_match, dc2_truth_match, mask = util.get_match(dc2_det, dc2_truth)
print(dc2_truth_match.dtype.names)
print(dc2_det_match.dtype.names)
print(len(dc2_truth_match))
print(len(dc2_det_match))
print(dc2_truth_match['ra'][0], dc2_truth_match['dec'][0])
print(dc2_det_match['alphawin_j2000'][0], dc2_det_match['deltawin_j2000'][0])
print(dc2_truth_match['ra'][1], dc2_truth_match['dec'][1])
print(dc2_det_match['alphawin_j2000'][1], dc2_det_match['deltawin_j2000'][1])