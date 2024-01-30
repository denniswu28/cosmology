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

# Read in all objects
dc2_truth = util.read_truth_catalog('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz')
dc2_truth = dc2_truth[dc2_truth['x']>0] # select valid
dc2_truth = dc2_truth[dc2_truth['gal_star'] == 0]
dc2_truth = dc2_truth[np.logical_and(np.logical_and(dc2_truth['ra'] < ra_max, dc2_truth['ra'] > ra_min),np.logical_and(dc2_truth['dec'] < dec_max, dc2_truth['dec'] > dec_min))]


for i in range(20,26):
    mag_threshold = i

    dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < mag_threshold]

    dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()
    dc2_thre_shear = np.copy(dc2_truth_shear[dc2_truth_thre['ind']])
    # dc2_truth_shear = dc2_truth_shear[np.logical_and(np.logical_and(dc2_truth_shear['ra'] < np.deg2rad(ra_max), dc2_truth_shear['ra'] > np.deg2rad(ra_min)),np.logical_and(dc2_truth_shear['dec'] < np.deg2rad(dec_max), dc2_truth_shear['dec'] > np.deg2rad(dec_min)))]
    
    dc2_truth_shear_limited = np.copy(dc2_truth_shear[dc2_truth['ind']])
    dc2_thre_shear = dc2_thre_shear[np.logical_and(np.logical_and(dc2_thre_shear['ra'] < np.deg2rad(ra_max), dc2_thre_shear['ra'] > np.deg2rad(ra_min)),np.logical_and(dc2_thre_shear['dec'] < np.deg2rad(dec_max), dc2_thre_shear['dec'] > np.deg2rad(dec_min)))]

    if i == 20:
        s1 = np.zeros(len(dc2_truth_shear_limited))
        s2 = np.zeros(len(dc2_truth_shear_limited))
        q = np.median(dc2_truth_shear_limited['q'],axis=1)
        beta = np.median(dc2_truth_shear_limited['pa'],axis=1)
        for i in range(len(dc2_truth_shear_limited)):
            try:
                s = galsim.Shear(q=1./q[i], beta=(90.+beta[i])*galsim.degrees)
                s = galsim._Shear(complex(s.g1,-s.g2))
                g1 = dc2_truth_shear_limited['g1'][i]/(1. - dc2_truth_shear_limited['k'][i])
                g2 = -dc2_truth_shear_limited['g2'][i]/(1. - dc2_truth_shear_limited['k'][i])
                g = galsim.Shear(g1=g1,g2=g2)
                s = g+s
                s1[i]=s.g1
                s2[i]=s.g2
            except:
                pass

        with open('data/s1.csv', 'x') as fs1:
            writer = csv.writer(fs1)
            writer.writerow(s1)
        del(s1)

        with open('data/s2.csv', 'x') as fs2:
            writer = csv.writer(fs2)
            writer.writerow(s2)
        del(s2)	

    s1_thre = np.zeros(len(dc2_thre_shear))
    s2_thre = np.zeros(len(dc2_thre_shear))
    q_thre = np.median(dc2_thre_shear['q'],axis=1)
    beta_thre = np.median(dc2_thre_shear['pa'],axis=1)
    for i in range(len(dc2_thre_shear)):
        try:
            s_thre = galsim.Shear(q=1./q_thre[i], beta=(90.+beta_thre[i])*galsim.degrees)
            s_thre = galsim._Shear(complex(s_thre.g1,-s_thre.g2))
            g1_thre = dc2_thre_shear['g1'][i]/(1. - dc2_thre_shear['k'][i])
            g2_thre = -dc2_thre_shear['g2'][i]/(1. - dc2_thre_shear['k'][i])
            g_thre = galsim.Shear(g1=g1_thre,g2=g2_thre)
            s_thre = g_thre+s_thre
            s1_thre[i]=s_thre.g1
            s2_thre[i]=s_thre.g2
        except:
            pass

    with open('data/s1_thre_'+ str(mag_threshold) +'.csv', 'x') as fs1_thre:
        writer = csv.writer(fs1_thre)
        writer.writerow(s1_thre)
        
    del(s1_thre)

    with open('data/s2_thre_'+ str(mag_threshold) +'.csv', 'x') as fs2_thre:
        writer = csv.writer(fs2_thre)
        writer.writerow(s2_thre)
        
    del(s2_thre)












