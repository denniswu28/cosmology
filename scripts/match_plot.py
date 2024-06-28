"""
This script is to compute the 2-point shear-shear correlation of the threshholded truth catalog using treecorr. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import numpy as np
import utilities as util
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt

fr = ['Y106','J129','H158','F184']

fr_num = 2
dec_min = -42.06
dec_max = -42
ra_min = 51
ra_max = 51.060
max_mag = 30
max_dec = -39
max_ra = 52

# Read in and cut truth objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
print(len(dc2_truth))
# for obj in dc2_truth:
#     if obj['mag_H158'] < max_mag:
#         max_mag = obj['mag_H158']
#         max_ra = obj['ra']
#         max_dec = obj['dec']

# dec_min = max_dec - 0.008
# dec_max = max_dec + 0.008
# ra_min = max_ra - 0.01
# ra_max = max_ra + 0.01

dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# Read in and cut dectection objects
dc2_det = util.read_det_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz')
dc2_det = dc2_det[dc2_det['number']>0] # select valid
dc2_det = dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5
dc2_det = dc2_det[(dc2_det['alphawin_j2000']>ra_min)&(dc2_det['alphawin_j2000']<ra_max)&(dc2_det['deltawin_j2000']>dec_min)&(dc2_det['deltawin_j2000']<dec_max)]

# Match both catalogs
dc2_det_match, dc2_truth_match = util.get_match(dc2_det, dc2_truth)

# Extract position data of matched detection catalog
ra_detection = dc2_det_match['alphawin_j2000']
dec_detection = dc2_det_match['deltawin_j2000']

ra_det = dc2_det['alphawin_j2000']
dec_det = dc2_det['deltawin_j2000']


# Extract position data of matched truth catalog
ra_truth = dc2_truth_match['ra']
dec_truth = dc2_truth_match['dec']

ra_tru = dc2_truth['ra']
dec_tru = dc2_truth['dec']



print(len(ra_truth))
print(len(ra_tru))

print(len(ra_det))

# plot RA/DEC dist of sources
plt.plot(ra_tru, dec_tru, linestyle="",marker=".", color='navy', label = r'Truth$')
plt.plot(ra_det, dec_det, linestyle="",marker=".", color='red', label = r'Det$', markersize = 1)
plt.xlabel(r'$\theta$ (deg)')
plt.ylabel(r'$\theta$ (deg)')
plt.legend()
plt.title("Position of matched detection and truth sources")
plt.savefig('sources_ori.pdf')
plt.clf()

# plot RA/DEC dist of sources
plt.plot(ra_tru, dec_tru, linestyle="",marker=".", color='navy', label = r'Truth$')
plt.plot(ra_detection, dec_detection, linestyle="",marker=".", color='red', label = r'Det$', markersize = 1)
plt.xlabel(r'$\theta$ (deg)')
plt.ylabel(r'$\theta$ (deg)')
plt.legend()
plt.title("Position of matched detection and truth sources")
plt.savefig('sources.pdf')
plt.clf()