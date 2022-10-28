"""
This script is to make a histogram of <H158 magnitude - object number density> from the truth catalog and Roman detection catalog.

Author: Dennis Wu
"""

# All the libs needed
import fitsio as fio
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


# Read in truth and detection catalog and get magnitudes

'''
truth_dir = '/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/truth/coadd'
det_dir = '/hpc/group/cosmology/phy-lsst/public/dc2_sim_output/detection'
true_mag = np.zeros([4, 1000000000])
det_mag = np.zeros([4, 1000000000])
fr = ['Y106','J129','H158','F184']

marker = 0
for f in np.sort(glob.glob(os.path.join(truth_dir, 'dc2_index_*.fits.gz'))):
	try:
		tmp = fio.FITS(f)[-1].read()
		for i in range(4):
			for j, mag_j in enumerate(tmp['mag_'+fr[i]]):
				true_mag[i][j+marker] = mag_j
		marker = marker + tmp['mag_'+fr[i]].size
	except:
		pass
	
	
marker = 0
for f in np.sort(glob.glob(os.path.join(det_dir, 'dc2_det_*.fits.gz'))):
	try:
		tmp = fio.FITS(f)[-1].read()
		for i in range(4):
			for j, mag_j in enumerate(tmp['mag_auto_'+fr[i]]):
				det_mag[i][j+marker] = mag_j
		marker = marker + tmp['mag_auto_'+fr[i]].size
	except:
		pass
'''


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
    except:
        print('-----fail'+f)
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


fr = ['Y106','J129','H158','F184']
# Plot the histogram
for i in range(4):
	n1, bins1, patches1 = plt.hist(rtd['mag_'+fr[i]], bins=np.linspace(10,35,1300), histtype='step', label='All Input')
	n2, bins2, patches2 = plt.hist(rd['mag_auto_'+fr[i]], bins=np.linspace(10,35,1300), histtype='step', label='Detected')
	plt.xlim( [10,35] )
	plt.yscale('log')
	plt.xlabel(fr[i]+' Magnitude')
	plt.ylabel(r'Number Density')
	plt.legend(loc='upper left')
	# Save both functions in a pdf file
	plt.savefig('hist_'+fr[i]+'.pdf')
	plt.clf()
	threshold = 0.9
	mag_threshold = 0.0
	for j in range(n1.size):
		if(n2[j] < n1[j]*threshold):
			mag_threshold = bins1[j]
	print('The threshold for magnitude is ', mag_threshold) #30.86220169361047 for 0.9 F184
