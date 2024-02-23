"""
This script is to compute the 2-point point-shear correlation of the threshholded truth catalog using treecorr. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
import pandas as pd
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
mag_threshold = 21

''' Write a dummy lens and source nofz to twopoint fits file. '''
# Read in and cut dectection objects
dc2_det = util.read_det_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz')
dc2_det = dc2_det[dc2_det['number']>0] # select valid
dc2_det = dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>18)] # select flag=0 and flux/err ratio >5
dc2_det = dc2_det[(dc2_det['alphawin_j2000']>ra_min)&(dc2_det['alphawin_j2000']<ra_max)&(dc2_det['deltawin_j2000']>dec_min)&(dc2_det['deltawin_j2000']<dec_max)] # boundary limit

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0] # select galaxy
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# Match both catalogs
dc2_det_match, dc2_truth_match = util.get_match(dc2_det, dc2_truth)

# Single out the lens
dc2_truth_thre = dc2_truth_match[dc2_truth_match['mag_'+fr[2]] < mag_threshold]
# dc2_truth_thre = dc2_truth_match

# Output the source, lens number
print("Matched source number: ", dc2_truth_match['ra'].size)
print("Matched lens number w/ m_H<" + str(mag_threshold) + " : ", dc2_truth_thre['ra'].size)

# Fetch nz dist from truth, create hist
ddir = '/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits'
nz_lens = util.nz_dist(ddir, dc2_truth_thre['ind'])
nz_source = util.nz_dist(ddir, dc2_truth_match['ind'])

hist_lens, bins_lens = np.histogram(nz_lens, bins = 299, range=(0.01,3.01), density = True)

hist_source, bins_source = np.histogram(nz_source, bins = 299, range=(0.01,3.01), density = True)

plt.stairs(hist_lens, bins_lens, label = 'Lens')
plt.stairs(hist_source, bins_source, label = 'Source')
plt.xlabel(r'redshift (z)')
plt.ylabel(r'Num Density')
plt.legend()
plt.title("Redshift Distribution of detected (matched) LSST galaxy")
plt.savefig('nz.pdf')
plt.clf()

bins_lens  = bins_lens[0:-1]
bins_source  = bins_source[0:-1]

# Format output file for CosmoCov
df = pd.DataFrame(np.transpose(np.array([bins_source, hist_source])), columns = ['z_min','n'])
df.to_csv('data/dc2_source_' + str(mag_threshold) + '.nz', sep = ' ', header = False, index = False, float_format = '{:.6e}'.format)
df = pd.DataFrame(np.transpose(np.array([bins_lens, hist_lens])), columns = ['z_min','n'])
df.to_csv('data/dc2_lens_' + str(mag_threshold) + '.nz', sep = ' ', header = False, index = False, float_format = '{:.6e}'.format)