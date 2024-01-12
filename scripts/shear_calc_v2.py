"""
This script is to calculate the reduced shear of both true catalog and selected bright objects within a certain threshhold. The lines reading in the catalog and calculating reduced shear is adapted from Prof. Troxel's script.

Author: Dennis Wu
"""

# All the libs needed
import fitsio as fio
import csv
import numpy as np
import galsim
import utilities as util

fr = ['Y106','J129','H158','F184']
dec_min = -42
dec_max = -38
ra_min = 51
ra_max = 56


# Read in all objects
dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()
dc2_truth_shear = dc2_truth_shear[np.logical_and(np.logical_and(dc2_truth_shear['ra'] < np.deg2rad(ra_max), dc2_truth_shear['ra'] > np.deg2rad(ra_min)),np.logical_and(dc2_truth_shear['dec'] < np.deg2rad(dec_max), dc2_truth_shear['dec'] > np.deg2rad(dec_min)))]

s1 = np.zeros(len(dc2_truth_shear))
s2 = np.zeros(len(dc2_truth_shear))
q = np.median(dc2_truth_shear['q'],axis=1)
beta = np.median(dc2_truth_shear['pa'],axis=1)

for i in range(len(dc2_truth_shear)):
    try:
        s = galsim.Shear(q=1./q[i], beta=(90.+beta[i])*galsim.degrees)
        s = galsim._Shear(complex(s.g1,-s.g2))
        g1 = dc2_truth_shear['g1'][i]/(1. - dc2_truth_shear['k'][i])
        g2 = -dc2_truth_shear['g2'][i]/(1. - dc2_truth_shear['k'][i])
        g = galsim.Shear(g1=g1,g2=g2)
        s = g+s
        s1[i]=s.g1
        s2[i]=s.g2
    except:
        pass

with open('data/s1_NG.csv', 'w') as fs1:
	writer = csv.writer(fs1)
	writer.writerow(s1)
	
del(s1)

with open('data/s2_NG.csv', 'w') as fs2:
	writer = csv.writer(fs2)
	writer.writerow(s2)
	
del(s2)
