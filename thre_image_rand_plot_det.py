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
dc2_det = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_det.fits.gz')[-1].read(columns=['mag_auto_J129','mag_auto_F184','mag_auto_H158','mag_auto_Y106', 'number', 'alphawin_j2000', 'deltawin_j2000','flux_auto','fluxerr_auto','flags'])
dc2_det=dc2_det[dc2_det['number']>0] # select positive
dc2_det=dc2_det[(dc2_det['flags']==0)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5

# Select objects based on threshold
dc2_det_thre = dc2_det[dc2_det['mag_auto_'+fr[2]] < 23]

# Extract position and shear data
ra_thre = dc2_det_thre['alphawin_j2000']
dec_thre = dc2_det_thre['deltawin_j2000']

del(dc2_det)
del(dc2_det_thre)

# Read the random catalog to arrays
with open('ra_rand_det.csv', newline='') as rafile:
    ra_rand = np.float_(list(csv.reader(rafile, delimiter=',')))

with open('dec_rand_det.csv', newline='') as decfile:
    dec_rand = np.float_(list(csv.reader(decfile, delimiter=',')))

ra_rand = ra_rand[0]
dec_rand = dec_rand[0]

# ind = np.random.randint(len(ra_rand), size = 1000)
ind = np.random.choice(len(ra_rand), size = 1000, replace=False)
ind_thre = np.random.choice(len(ra_thre), size = 1000, replace=False)

ra_rand = ra_rand[ind]
dec_rand = dec_rand[ind]
ra_thre = ra_thre[ind_thre]
dec_thre = dec_thre[ind_thre]


plt.plot(ra_rand, dec_rand, linestyle="",marker="*", color='blue', label = r'det rand catalog objects')
plt.xlabel(r'Ascension (degree)')
plt.ylabel(r'Declination (degree)')
plt.legend()
plt.title("ra-dec coordinates plot for rand catalog")

plt.savefig('rand_coordinates_det.pdf')
plt.clf()

plt.plot(ra_thre, dec_thre, linestyle="",marker="*", color='red', label = r'det thre catalog objects')
plt.xlabel(r'Ascension (degree)')
plt.ylabel(r'Declination (degree)')
plt.legend()
plt.title("ra-dec coordinates plot for thre catalog")

plt.savefig('thre_coordinates_det.pdf')
plt.clf()


