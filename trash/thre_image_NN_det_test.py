"""
This script is to compute the 2-point count-count correlation of the threshholded detection catalog using treecorr with Landy-Szalay estimator. The lines reading in the catalog is adapted from Prof. Troxel's script.

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
import glob


fr = ['Y106','J129','H158','F184']
mag_threshold = 23

# Read in all objects
# dc2_det = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_det.fits.gz')[-1].read(columns=['mag_auto_J129','mag_auto_F184','mag_auto_H158','mag_auto_Y106', 'number', 'alphawin_j2000', 'deltawin_j2000','flux_auto','fluxerr_auto','flags'])

start  = 0
dc2_det = None
for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz'))[5:-4]):
	try:
		tmp = fio.FITS(f)[-1].read(columns=['mag_auto_J129','mag_auto_F184','mag_auto_H158','mag_auto_Y106', 'number', 'alphawin_j2000', 'deltawin_j2000','flux_auto','fluxerr_auto','flags'])
		if dc2_det is None:
			dc2_det = np.zeros(100000000,dtype=tmp.dtype)
		for col in dc2_det.dtype.names:
			dc2_det[col][start:start+len(tmp)] = tmp[col]
		start+=len(tmp)
	except:
		print('-----fail'+f)
		pass

dc2_det=dc2_det[dc2_det['number']>0] # select positive
dc2_det=dc2_det[(dc2_det['flags']<4)&(dc2_det['flux_auto']/dc2_det['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5

# Select objects based on threshold
# dc2_det_thre = dc2_det[dc2_det['mag_auto_'+fr[2]] < mag_threshold]

# Extract position and shear data
ra_thre = dc2_det['alphawin_j2000']
dec_thre = dc2_det['deltawin_j2000']

# Release RAM
del(dc2_det)


# Read the random catalog to arrays
# with open('ra_rand_det.csv', newline='') as rafile:
#     ra_rand = np.float_(list(csv.reader(rafile, delimiter=',')))
# ra_rand = ra_rand[0]

# with open('dec_rand_det.csv', newline='') as decfile:
#     dec_rand = np.float_(list(csv.reader(decfile, delimiter=',')))
# dec_rand = dec_rand[0]

# ra_rand = np.deg2rad(ra_rand)
# dec_rand = np.deg2rad(dec_rand)
# ra_thre = np.deg2rad(ra_thre)
# dec_thre = np.deg2rad(dec_thre)

# Generate the treecorr catalogs
N = 10
cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, ra_units='degrees', dec_units='degrees')
# cat_rand = treecorr.Catalog(ra = ra_rand, dec = dec_rand, ra_units='degrees', dec_units='degrees', patch_centers=cat_thre.patch_centers)

# Construct the count-shear corr object
dd = treecorr.NNCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin')

# Compute the count-count corr
dd.process(cat_thre)
# dr.process(cat_thre, cat_rand)
# rr.process(cat_rand)
# dd.calculateXi()
# dd.write("nn_corr_det.fits")

# plot the corr functions for NN
# r = np.exp(dd.meanlogr)
# xi = dd.xi
# sig = np.sqrt(dd.varxi)

print(dd.npairs)

# plt.errorbar(r, xi, yerr=sig, linestyle="",marker="*", color='blue', label = r'$\xi_(\theta)$')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$\theta$ (arcmin)')
# plt.ylabel(r'$\xi$')
# plt.legend()
# plt.title("2-point NN correlation of det catalog")

# plt.savefig('corr_func_threshold_NN_det.pdf')
# plt.clf()

