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



fr = ['Y106','J129','H158','F184']

# Read in all objects

start  = 0
rtd = None
for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd/dc2_index_*.fits.gz'))):
    print(i)
    try:
        tmp = fio.FITS(f)[-1].read()
        if rtd is None:
            print('test')
            rtd = np.zeros(100000000,dtype=tmp.dtype)
        for col in rtd.dtype.names:
            rtd[col][start:start+len(tmp)] = tmp[col]
        start+=len(tmp)
        print('-----success: '+f)
    except:
        print('-----fail: '+f)
        pass
        

start  = 0
rd = None
for i,f in enumerate(np.sort(glob.glob('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection/dc2_det_*.fits.gz'))):
    # if i in mask:
    tmp = fio.FITS(f)[-1].read()
    if rd is None:
        rd = np.zeros(100000000,dtype=tmp.dtype)
    for col in rd.dtype.names:
        rd[col][start:start+len(tmp)] = tmp[col]
    # mask.append(i)
    start+=len(tmp)


rd=rd[rd['number']>0] # select positive
rd=rd[(rd['flags']==0)&(rd['flux_auto']/rd['fluxerr_auto']>5)] # select flag=0 and flux/err ratio >5


# calculate threshold

threshold = 0.85
mag_threshold = np.zeros(4)
# Plot the histogram
for i in range(4):
	n1, bins1, patches1 = plt.hist(rtd['mag_'+fr[i]], bins=np.linspace(10,35,1300), histtype='step', label='All Input')
	n2, bins2, patches2 = plt.hist(rd['mag_auto_'+fr[i]], bins=np.linspace(10,35,1300), histtype='step', label='Detected')
	for j in range(n1.size):
		if(n2[j] > n1[j]*threshold and n2[j] > 1000):
			mag_threshold[i] = bins1[j]
	print('The threshold at ' + fr[i] + ' is ', mag_threshold[i])
	plt.clf()
	


# Select objects based on threshold

dc2_truth = fio.FITS('/hpc/group/cosmology/phy-lsst/public/andy_tru.fits.gz')[-1].read()

#for i in range(4):
	#dc2_truth = dc2_truth[dc2_truth['mag_'+fr[i]] < mag_threshold[i]]
	
dc2_truth = dc2_truth[dc2_truth['mag_'+fr[2]] < mag_threshold[2]]

dc2_truth_shear = fio.FITS('/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/dc2_truth_gal.fits')[-1].read()

dc2_thre_shear = np.copy(dc2_truth_shear[dc2_truth['ind']])

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

with open('s1.csv', 'x') as fs1:
	writer = csv.writer(fs1)
	writer.writerow(s1)
	
del(s1)

with open('s2.csv', 'x') as fs2:
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



with open('s1_thre.csv', 'x') as fs1_thre:
	writer = csv.writer(fs1_thre)
	writer.writerow(s1_thre)
	
del(s1_thre)

with open('s2_thre.csv', 'x') as fs2_thre:
	writer = csv.writer(fs2_thre)
	writer.writerow(s2_thre)
	
del(s2_thre)












