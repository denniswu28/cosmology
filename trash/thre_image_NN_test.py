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
import treecorr
import glob


fr = ['Y106','J129','H158','F184']
mag_threshold = 23

# Read in all objects
# dc2_truth = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_tru.fits.gz')[-1].read(columns=['mag_J129','mag_F184','mag_H158','mag_Y106', 'ind', 'ra', 'dec'])

start  = 0
dc2_truth = None
for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz'))):
    # print(i)
    try:
        tmp = fio.FITS(f)[-1].read(columns=['ra','dec', 'mag_J129','mag_F184','mag_H158','mag_Y106', 'ind'])
        if dc2_truth is None:
            print('test')
            dc2_truth = np.zeros(100000000,dtype=tmp.dtype)
        for col in dc2_truth.dtype.names:
            dc2_truth[col][start:start+len(tmp)] = tmp[col]
        start+=len(tmp)
    except:
        print('-----fail'+f)
        pass

# Select objects based on threshold
# dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < mag_threshold]

# # Read in shear data
# dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read(columns=['ra','dec'])

# Create new array with only data of bright objects 
# dc2_thre_shear = np.copy(dc2_truth_shear[dc2_truth_thre['ind']])

# Limit the objects in the survey area
# dec_min = -0.6691307
# dec_max = -0.6156614
# ra_min = 0.890117
# ra_max = 0.977386

dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]

# Extract position and shear data
ra_thre = dc2_truth['ra']
dec_thre = dc2_truth['dec']

# Release RAM
del(dc2_truth)
# del(dc2_thre_shear)


# # Read the random catalog to arrays
# with open('ra_rand_man.csv', newline='') as rafile:
#     ra_rand = np.float_(list(csv.reader(rafile, delimiter=',')))
# ra_rand = ra_rand[0]

# with open('dec_rand_man.csv', newline='') as decfile:
#     dec_rand = np.float_(list(csv.reader(decfile, delimiter=',')))
# dec_rand = dec_rand[0]

# Read the random catalog to arrays
# with open('ra_rand_man_1.csv', newline='') as rafile:
#     ra_rand = np.float_(list(csv.reader(rafile, delimiter=',')))
# ra_rand = ra_rand[0]

# with open('dec_rand_man_1.csv', newline='') as decfile:
#     dec_rand = np.float_(list(csv.reader(decfile, delimiter=',')))
# dec_rand = dec_rand[0]

# Generate the treecorr catalogs and patches with k-means
cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, ra_units='degrees', dec_units='degrees')
#cat_rand = treecorr.Catalog(ra = ra_rand, dec = dec_rand, ra_units='degrees', dec_units='degrees', patch_centers=cat_thre.patch_centers)

# Construct the count-shear corr object
dd = treecorr.NNCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin')
# dr = treecorr.NNCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin')
# rr = treecorr.NNCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin')

# Compute the count-count corr
dd.process(cat_thre)
# dr.process(cat_thre, cat_rand)
# rr.process(cat_rand)
# dd.calculateXi()
# dd.write("nn_corr_1.fits")

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
# plt.title("2-point NN correlation")

# plt.savefig('corr_func_threshold_NN.pdf')
# plt.clf()





