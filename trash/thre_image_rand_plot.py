"""
This script is to compute the 2-point count-count correlation of the threshholded truth catalog using treecorr with Landy-Szalay estimator. The lines reading in the catalog is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
import matplotlib
matplotlib.use ('agg')
import matplotlib.pyplot as plt
import healpy as hp



fr = ['Y106','J129','H158','F184']

# Read in all objects
dc2_truth = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_tru.fits.gz')[-1].read(columns=['mag_J129','mag_F184','mag_H158','mag_Y106', 'ind'])

# Select objects based on threshold
dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < 23]

# Read in shear data
dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read(columns=['ra','dec'])

# Create new array with only data of bright objects
dc2_thre_shear = np.copy(dc2_truth_shear[dc2_truth_thre['ind']])

# Limit the objects in the survey area
dec_min = -0.6691307
dec_max = -0.6156614
ra_min = 0.890117
ra_max = 0.977386
dc2_thre_shear = dc2_thre_shear[np.logical_and(np.logical_and(dc2_thre_shear['ra'] < ra_max, dc2_thre_shear['ra'] > ra_min), \
                                np.logical_and(dc2_thre_shear['dec'] < dec_max, dc2_thre_shear['dec'] > dec_min))]


# Extract position and shear data
ra_thre = dc2_thre_shear['ra']
dec_thre = dc2_thre_shear['dec']

del(dc2_truth)
del(dc2_thre_shear)
del(dc2_truth_thre)

# Read the random catalog to arrays
with open('ra_rand_man.csv', newline='') as rafile:
    ra_rand = np.float_(list(csv.reader(rafile, delimiter=',')))

with open('dec_rand_man.csv', newline='') as decfile:
    dec_rand = np.float_(list(csv.reader(decfile, delimiter=',')))

ra_rand = ra_rand[0]
dec_rand = dec_rand[0]

ind = np.random.choice(len(ra_rand), size = 1000, replace=False)
ind_thre = np.random.choice(len(ra_thre), size = 1000, replace=False)

ra_rand = np.rad2deg(ra_rand[ind])
dec_rand = np.rad2deg(dec_rand[ind])
ra_thre = np.rad2deg(ra_thre[ind_thre])
dec_thre = np.rad2deg(dec_thre[ind_thre])

plt.plot(ra_rand, dec_rand, linestyle="",marker="*", color='blue', label = r'truth rand catalog objects')
plt.xlabel(r'Ascension (degree)')
plt.ylabel(r'Declination (degree)')
plt.legend()
plt.title("ra-dec coordinates plot for rand catalog")

plt.savefig('rand_coordinates.pdf')
plt.clf()

plt.plot(ra_thre, dec_thre, linestyle="",marker="*", color='red', label = r'truth thre catalog objects')
plt.xlabel(r'Ascension (degree)')
plt.ylabel(r'Declination (degree)')
plt.legend()
plt.title("ra-dec coordinates plot for thre catalog")

plt.savefig('thre_coordinates.pdf')
plt.clf()


