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

with open('s1.csv') as fs1:
	reader = csv.reader(fs1, delimiter=',')
	for row in reader:
		s1 = row


with open('s2.csv') as fs2:
	reader = csv.reader(fs2, delimiter=',')
	for row in reader:
		s2 = row


ra = dc2_truth_shear['ra']
dec = dc2_truth_shear['dec']
#g1 = dc2_truth_shear['g1']
#g2 = -dc2_truth_shear['g2']

s1_thre = np.zeros(len(dc2_thre_shear))
s2_thre = np.zeros(len(dc2_thre_shear))


with open('s1_thre.csv') as fs1_thre:
	reader = csv.reader(fs1_thre, delimiter=',')
	for row in reader:
		s1_thre = row

	
with open('s2_thre.csv') as fs2_thre:
	reader = csv.reader(fs2_thre, delimiter=',')
	for row in reader:
		s2_thre = row

ra_thre = dc2_thre_shear['ra']
dec_thre = dc2_thre_shear['dec']
#g1_thre = dc2_thre_shear['g1']
#g2_thre = -dc2_thre_shear['g2']

print("total number: ", ra.size)
print("thresholded number: ", ra_thre.size)

# Generate the treecorr catalog
cat = treecorr.Catalog(ra = ra, dec = dec, g1 = s1, g2 = s2, ra_units='radians', dec_units='radians')
cat_thre = treecorr.Catalog(ra = ra_thre, dec = dec_thre, g1 = s1_thre, g2 = s2_thre, ra_units='radians', dec_units='radians')

# Construct the shear-shear corr object
gg = treecorr.GGCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin')
gg_thre = treecorr.GGCorrelation(min_sep=2.5, max_sep=250, nbins=20, sep_units='arcmin')
# Compute the shear-shear corr
gg.process(cat)
gg_thre.process(cat_thre)

# Plot the xi functions
r = np.exp(gg.meanlogr)
xip = gg.xip
xim = gg.xim
sig = np.sqrt(gg.varxip)

plt.plot(r, xip, color='blue', linestyle="",marker=".", label = r'$\xi_+(\theta)$')
lp = plt.errorbar(-r, xip, yerr=sig, color='blue')

plt.plot(r, xim, color='green', linestyle="",marker=".", label = r'$\xi_-(\theta)$')
lm = plt.errorbar(-r, xim, yerr=sig, color='green')

plt.errorbar(r[xip>0], xip[xip>0], yerr=sig[xip>0], color='k', lw=0.1, ls='')
plt.errorbar(r[xip<0], -xip[xip<0], yerr=sig[xip<0], color='k', lw=0.1, ls='')


r_thre = np.exp(gg_thre.meanlogr)
xip_thre = gg_thre.xip
xim_thre = gg_thre.xim
sig_thre = np.sqrt(gg_thre.varxip)

plt.plot(r_thre, xip_thre, color='blue', linestyle="",marker="*", label = r'$\xi_{+}^{\prime} (\theta)$')
lp_thre = plt.errorbar(-r_thre, xip_thre, yerr=sig_thre, color='blue')

plt.plot(r_thre, xim_thre, color='green', linestyle="",marker="*", label = r'$\xi_{-}^{\prime} (\theta)$')
lm_thre = plt.errorbar(-r_thre, xim_thre, yerr=sig_thre, color='green')

plt.errorbar(r_thre[xip_thre>0], xip_thre[xip_thre>0], yerr=sig_thre[xip_thre>0], color='k', lw=0.1, ls='')
plt.errorbar(r_thre[xip_thre<0], -xip_thre[xip_thre<0], yerr=sig_thre[xip_thre<0], color='k', lw=0.1, ls='')

plt.grid(color = 'k', linestyle = '-', linewidth = 1)

plt.xscale('log')
plt.yscale('log')
#plt.xlim([1,500])
plt.xlabel(r'$\theta$ (arcmin)')

plt.legend()
plt.ylabel(r'$\xi_{+,-}$')

# Save both functions in a pdf file
plt.savefig('corr_func_threshold.pdf')
plt.clf()




















