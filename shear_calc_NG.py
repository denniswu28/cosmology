"""
This script is to calculate the reduced shear of both true catalog and selected bright objects within a certain threshhold. The lines reading in the catalog and calculating reduced shear is adapted from Prof. Troxel's script.

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



fr = ['Y106','J129','H158','F184']


# Read in all objects
dc2_truth = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_tru.fits.gz')[-1].read()

# Select objects based on threshold
dc2_truth_thre = dc2_truth[dc2_truth['mag_'+fr[2]] < 24]

dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()

dc2_thre_shear = np.copy(dc2_truth_shear[dc2_truth_thre['ind']])

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

with open('s1_NG.csv', 'w') as fs1:
	writer = csv.writer(fs1)
	writer.writerow(s1)
	
del(s1)

with open('s2_NG.csv', 'w') as fs2:
	writer = csv.writer(fs2)
	writer.writerow(s2)
	
del(s2)	


# s1_thre = np.zeros(len(dc2_thre_shear))
# s2_thre = np.zeros(len(dc2_thre_shear))
# q_thre = np.median(dc2_thre_shear['q'],axis=1)
# beta_thre = np.median(dc2_thre_shear['pa'],axis=1)
# for i in range(len(dc2_thre_shear)):
#     try:
#         s_thre = galsim.Shear(q=1./q_thre[i], beta=(90.+beta_thre[i])*galsim.degrees)
#         s_thre = galsim._Shear(complex(s_thre.g1,-s_thre.g2))
#         g1_thre = dc2_thre_shear['g1'][i]/(1. - dc2_thre_shear['k'][i])
#         g2_thre = -dc2_thre_shear['g2'][i]/(1. - dc2_thre_shear['k'][i])
#         g_thre = galsim.Shear(g1=g1_thre,g2=g2_thre)
#         s_thre = g_thre+s_thre
#         s1_thre[i]=s_thre.g1
#         s2_thre[i]=s_thre.g2
#     except:
#         pass



# with open('s1_thre_NG.csv', 'w') as fs1_thre:
# 	writer = csv.writer(fs1_thre)
# 	writer.writerow(s1_thre)
	
# del(s1_thre)

# with open('s2_thre_NG.csv', 'w') as fs2_thre:
# 	writer = csv.writer(fs2_thre)
# 	writer.writerow(s2_thre)
	
# del(s2_thre)












